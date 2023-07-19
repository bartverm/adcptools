classdef TidalModel < DataModel
    % Velocity model to include tidal dynamics
    %
    %   TidalModel properties:
    %   constituents
    %
    %   TidalModel methods:
    %   get_tidal_pars - compute the tidal amplitude and phase
    %
    %   see also: DataModel, TaylorModel

    properties
        % TidalModel/constituents constituents to fit data with
        %
        %   MxN array defining tidal constituents to be included in the model
        %   for data to be fitted. Every row is a different component, e.g. for
        %   velocity: row 1 is x component, row 2 is y component and row 3 is z
        %   component. For a scalar only one row is provided.
        %   Columns represent different constituents. Zero or NaN values are
        %   skipped. For every given value, two model parameters are fitted, which
        %   are the coefficients for the sin and cos function. From those amplitude
        %   and phases are computed. A residual component is always included.
        %
        %   see also: TidalModel, get_tidal_pars
        constituents
    end

    methods
        function M = get_model(obj, d_time, ~, ~, ~, ~)
            % This model fits the following parameters to the velocity
            % within each cell:
            % u = u_0 + sum_n (a_n cos(2pi/T_n * t) + b_n sin(2pi/T_n * t))
            % where n loops over all entered constituents (subtidal
            % constituent is always present)
            % Input:
            % d_time IN HOURS
            % Output:
            % Model matrices such that u = Mu*pars (roughly)
            npars = obj.npars;
            ncomp = obj.ncomponents;
            d_secs = seconds(d_time);

            max_pars = max(npars);
            M = nan(numel(d_time), max_pars, ncomp);
            M(:,1,:) = 1; %residual
            for c_comp = 1:ncomp
                n_const = sum(isfinite(obj.constituents(c_comp,:)) &...
                    obj.constituents(c_comp,:)~=0);
                for c_const = 1:n_const
                    M(:,2*c_const,c_comp) = ...
                        cos( ...
                        2*pi/obj.constituents(c_comp,c_const) *...
                        hours(d_secs) ...
                        );
                    M(:,2*c_const + 1,c_comp) = sin( ...
                        2*pi/obj.constituents(c_comp,c_const) *...
                        hours(d_secs) ...
                        );
                end
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
            %   the Solver class.
            %
            %   cov_pars holds the covariance matrix of the model parameters.
            %
            %   see also: TidalModel, Solver

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
           nconst = sum(isfinite(obj.constituents) &...
               obj.constituents~=0,2)';
           val = nconst*2 + 1;
        end
    end

end
