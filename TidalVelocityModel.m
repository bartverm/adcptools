classdef TidalVelocityModel < VelocityModel
    % Velocity model to include tidal dynamics
    %
    %   TidalVelocityModel properties:
    %   constituentsU - tidal constituents for x-component of velocity
    %   constituentsV - tidal constituents for y-component of velocity
    %   constituentsW - tidal constituents for z-component of velocity
    %
    %   TidalVelocityModel methods:
    %   get_tidal_pars - compute the tidal amplitude and phase
    %
    %   see also: VelocityModel, TaylorExpandedVelocity

    properties
        % TidalVelocityModel/constituentsU constituents for x-component of velocity
        %
        %   1xN row vector defining tidal constituents to be included in the model
        %   for the x-component of the velocity. Every string given indicates the
        %   name of the constituent to be included. For every given value, two
        %   model parameters are fitted, which are the coefficients for the cos and
        %   sin functions. From those amplitude and phases are computed. A residual
        %   is always included.
        %
        %   see also: TidalVelocityModel, constituentsV, constituentsW,
        %               get_tidal_pars
        constituentsU(1,:) cell = {};

        % TidalVelocityModel/constituentsV constituents for y-component of velocity
        %
        %   1xN row vector defining tidal constituents to be included in the model
        %   for the y-component of the velocity. Every string given indicates the
        %   name of the constituent to be included. For every given value, two
        %   model parameters are fitted, which are the coefficients for the cos and
        %   sin functions. From those amplitude and phases are computed. A residual
        %   is always included.
        %
        %   see also: TidalVelocityModel, constituentsU, constituentsW,
        %               get_tidal_pars
        constituentsV(1,:) cell = {};

        % TidalVelocityModel/constituentsW constituents for z-component of velocity
        %
        %   1xN row vector defining tidal constituents to be included in the model
        %   for the z-component of the velocity. Every string given indicates the
        %   name of the constituent to be included. For every given value, two
        %   model parameters are fitted, which are the coefficients for the cos and
        %   sin functions. From those amplitude and phases are computed. A residual
        %   is always included.
        %
        %   see also: TidalVelocityModel, constituentsU, constituentsV,
        %               get_tidal_pars
        constituentsW(1,:) cell = {};

    end

    methods
        function [Mu, Mv, Mw] = get_model(obj, d_time, d_s, d_n, d_z, d_sigma)
            % This model fits the following parameters to the velocity
            % within each cell:
            % u = u_0 + sum_n (a_n cos(2pi/T_n * t) + b_n sin(2pi/T_n * t))
            % where n loops over all entered constituents (subtidal
            % constituent is always present)
            % Input:
            % d_time IN HOURS
            % Output:
            % Model matrices such that u = Mu*pars (roughly)
            npars = [2*numel(obj.constituentsU) + 1, ...
                2*numel(obj.constituentsV) + 1, ...
                2*numel(obj.constituentsW) + 1];
            d_secs = seconds(d_time);
            Tu = obj.const_to_periods(obj.constituentsU);
            Tv = obj.const_to_periods(obj.constituentsV);
            Tw = obj.const_to_periods(obj.constituentsW);

            Mu = ones(numel(d_time), npars(1));
            for c = 1:numel(obj.constituentsU)
                Mu(:,2*c) = cos(2*pi/Tu(c) * hours(d_secs));
                Mu(:,2*c + 1) = sin(2*pi/Tu(c) * hours(d_secs));
            end
            Mv = ones(numel(d_time), npars(2));

            for c = 1:numel(obj.constituentsV)
                Mv(:,2*c) = cos(2*pi/Tv(c) * hours(d_secs));
                Mv(:,2*c + 1) = sin(2*pi/Tv(c) * hours(d_secs));
            end
            Mw = ones(numel(d_time), npars(3));

            for c = 1:numel(obj.constituentsW)
                Mw(:,2*c) = cos(2*pi/Tw(c) * hours(d_secs));
                Mw(:,2*c + 1) = sin(2*pi/Tw(c) * hours(d_secs));
            end

        end

        function [pars_h, cov_pars_h] = get_tidal_pars(obj, pars, cov_pars)
            % Compute tidal amplitude and phase from sin and cos coefficients
            %
            %   [pars_h, cov_pars_h] = get_tidal_pars(obj, pars, cov_pars) computes the
            %   amplitude and phases based on the sin and cos coefficients:
            %   u = A cos(omega t - phi) = a cos(omega t) + b sin(omega t) with:
            %   A = sqrt(a^2 + b^2)
            %   phi = arctan(b / a)
            %
            %   pars is an array holding the tidal model parameters as produced with
            %   the VelocitySolver class.
            %
            %   cov_pars holds the covariance matrix of the model parameters.
            %
            %   see also: TidalVelocityModel, VelocitySolver

            npars = obj.npars;
            subtidal_idx = [1, npars(1) + 1, sum(npars(1:2)) + 1];
            pars_h = zeros(size(pars));
            jac = zeros(size(cov_pars)); % jacobian of the transformation (subtidal independent of other under this transformations)
            idx = 1;
            while idx <= size(pars_h,2)
                if ismember(idx, subtidal_idx)
                    pars_h(:,idx) = pars(:,idx);
                    jac(:,idx,idx) = 1;
                    idx = idx + 1;
                else
                    pars_h(:,idx) = sqrt(pars(:,idx).^2 + pars(:,idx+1).^2); % Amplitude of constituent
                    pars_h(:,idx+1) = atan(pars(:,idx+1)./pars(:,idx));     % Phase of constituent
                    jac(:,idx,idx) = pars(:,idx)./pars_h(:,idx);            % d(sqrt(a^2 + b^2))/da = a/A
                    jac(:,idx, idx+1) = pars(:,idx+1)./pars_h(:,idx);       % d(sqrt(a^2 + b^2))/db = b/A
                    jac(:,idx+1,idx) = -pars(:,idx+1)./pars_h(:,idx).^2;    % d(atan(b/a))/da = -b/A^2
                    jac(:,idx+1,idx+1) = pars(:,idx)./pars_h(:,idx).^2;     % d(atan(b/a))/db = a/A^2
                    idx = idx + 2;
                end
            end
            cov_pars_h = helpers.matmult(cov_pars,permute(jac,[1,3,2]),2,3);   % cov_pars * J'
            cov_pars_h = helpers.matmult(jac,cov_pars_h,2,3);                  % J * (cov_pars * J')
        end
    end

    

    methods(Access=protected)
        function val = get_npars(obj)
            val = [2*numel(obj.constituentsU) + 1, ...
                2*numel(obj.constituentsV) + 1, ...
                2*numel(obj.constituentsW) + 1];
        end
    end

end
