classdef TaylorExpandedVelocity < VelocityModel
% Velocity model based on Taylor expansions.
%
%   TaylorExpandedVelocity properties:
%   s_order - expand to given orders in s coordinate
%   n_order - expand to given orders in n coordinate
%   z_order - expand to given orders in z coordinate
%   sigma_order - expand to given orders in sigma coordinate
%   time_order - expand to given orders in time
%
%   TaylorExpandedVelocity methods:
%   get_velocity - return velocity based on model parameters
%
%  see also: VelocityModel

    properties
        % TaylorExpandedVelocity/s-order order of expansion in s-direction
        %
        %   1x3 row vector holding non-negative integers that specify up
        %   to which order the velocity should be expanded in the
        %   s-direction for each of the three cartesian velocity components
        %   [2 1 0] means expand x-velocity up to second order in
        %   s-direction, y-velocity up to first order and z-velocity up to
        %   0th order. For the x-velocity this results in three model
        %   parameters (mean x-vel, dx-vel/ds, d^2 x_vel/ds^2), for the
        %   y-velocity in two parameters (mean y-vel, dy-vel/dx) and for
        %   z-velocity in 1 parameter (mean z-vel).
        %
        %   Default is: [0, 0, 0]
        %
        %   see also: TaylorExpandedVelocity, n_order, z_order,
        %   sigma_order, time_order
        s_order (1,3) double {...
            mustBeFinite, mustBeInteger, mustBeNonnegative} = [0 0 0];

        % TaylorExpandedVelocity/n-order order of expansion in n-direction
        %
        %   1x3 row vector holding non-negative integers that specify up
        %   to which order the velocity should be expanded in the
        %   n-direction for each of the three cartesian velocity
        %   components. 
        %
        %   Default is: [0, 0, 0]
        %
        %   See for an example: s_order
        %
        %   see also: TaylorExpandedVelocity, s_order, z_order,
        %   sigma_order, time_order
        n_order (1,3) double {...
            mustBeFinite, mustBeInteger, mustBeNonnegative} = [0 0 0];

        % TaylorExpandedVelocity/z-order order of expansion in z-direction
        %
        %   1x3 row vector holding non-negative integers that specify up
        %   to which order the velocity should be expanded in the
        %   z-direction for each of the three cartesian velocity
        %   components. 
        %
        %   Default is: [0, 0, 0]
        %
        %   See for an example: s_order
        %
        %   see also: TaylorExpandedVelocity, s_order, n_order,
        %   sigma_order, time_order        
        z_order (1,3) double {...
            mustBeFinite, mustBeInteger, mustBeNonnegative} = [0 0 0];

        % TaylorExpandedVelocity/sigma-order order of expansion in sigma direction
        %
        %   1x3 row vector holding non-negative integers that specify up
        %   to which order the velocity should be expanded in the
        %   sigma-direction for each of the three cartesian velocity
        %   components. 
        %
        %   Default is: [0, 0, 0]
        %
        %   See for an example: s_order
        %
        %   see also: TaylorExpandedVelocity, s_order, z_order,
        %   n_order, time_order
        sigma_order (1,3) double {...
            mustBeFinite, mustBeInteger, mustBeNonnegative} = [0 0 0];

        % TaylorExpandedVelocity/time-order order of expansion in time
        %
        %   1x3 row vector holding non-negative integers that specify up
        %   to which order the velocity should be expanded in 
        %   time for each of the three cartesian velocity
        %   components. 
        %
        %   Default is: [0, 0, 0]
        %
        %   See for an example: s_order
        %
        %   see also: TaylorExpandedVelocity, s_order, z_order,
        %   n_order, n_order
        time_order (1,3) double {...
            mustBeFinite, mustBeInteger, mustBeNonnegative} = [0 0 0];

    end
    methods(Access = protected)
        function val = get_npars(obj)
            val = ones(1,3);
            val = val +...
                obj.s_order +...
                obj.n_order +...
                obj.z_order +...
                obj.sigma_order +...
                obj.time_order;
        end
    end
    methods
        function [Mu, Mv, Mw]=get_model(obj, time, d_s, d_n, d_z, d_sigma)
            Mu = obj.combine_expands(1, time, d_s, d_n, d_z, d_sigma);
            Mv = obj.combine_expands(2, time, d_s, d_n, d_z, d_sigma);
            Mw = obj.combine_expands(3, time, d_s, d_n, d_z, d_sigma);
        end
        function M = combine_expands(obj, dim, time, d_s, d_n, d_z, d_sigma)
            M= [ones(numel(time),1),...
                obj.taylor_expand(time,obj.time_order(dim)),...
                obj.taylor_expand(d_s,obj.s_order(dim)),...
                obj.taylor_expand(d_n,obj.n_order(dim)),...
                obj.taylor_expand(d_z,obj.z_order(dim)),...
                obj.taylor_expand(d_sigma,obj.sigma_order(dim))];
        end
    end
    methods (Static)
        function M = taylor_expand(var, ord)
            nin = numel(var);
            M = cumprod(repmat(reshape(var, [], 1), 1, ord), 2)./... % (x-x0)^n
                cumprod( cumsum( ones( nin, ord), 2), 2); % n!
        end
    end
end


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
        %   for the x-component of the velocity. Every value given indicates the
        %   period of the constituent to be included. For every given value, two
        %   model parameters are fitted, which are the coefficients for the sin and
        %   cos function. From those amplitude and phases are computed. A residual
        %   is always included.
        %
        %   see also: TidalVelocityModel, constituentsV, constituentsW,
        %               get_tidal_pars
        constituentsU(1,:) double {mustBeFinite, mustBePositive} = []

        % TidalVelocityModel/constituentsV constituents for y-component of velocity
        %
        %   1xN row vector defining tidal constituents to be included in the model
        %   for the y-component of the velocity. Every value given indicates the
        %   period of the constituent to be included. For every given value, two
        %   model parameters are fitted, which are the coefficients for the sin and
        %   cos function. From those amplitude and phases are computed. A residual
        %   is always included.
        %
        %   see also: TidalVelocityModel, constituentsU, constituentsW,
        %               get_tidal_pars
        constituentsV(1,:) double {mustBeFinite, mustBePositive} = []

        % TidalVelocityModel/constituentsW constituents for z-component of velocity
        %
        %   1xN row vector defining tidal constituents to be included in the model
        %   for the z-component of the velocity. Every value given indicates the
        %   period of the constituent to be included. For every given value, two
        %   model parameters are fitted, which are the coefficients for the sin and
        %   cos function. From those amplitude and phases are computed. A residual
        %   is always included.
        %
        %   see also: TidalVelocityModel, constituentsU, constituentsV,
        %               get_tidal_pars
        constituentsW(1,:) double {mustBeFinite, mustBePositive} = []
    end

    methods
        function [Mu, Mv, Mw] = get_model(obj, d_time, ~, ~, ~, ~) %What about spatial variation?
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
            d_secs = seconds(d_time);

            Mu = ones(numel(d_time), npars(1));
            for c = 1:numel(obj.constituentsU)
                Mu(:,2*c) = cos(2*pi/obj.constituentsU(c) * hours(d_secs));
                Mu(:,2*c + 1) = sin(2*pi/obj.constituentsU(c) * hours(d_secs));
            end
            Mv = ones(numel(d_time), npars(2));

            for c = 1:numel(obj.constituentsV)
                Mv(:,2*c) = cos(2*pi/obj.constituentsV(c) * hours(d_secs));
                Mv(:,2*c + 1) = sin(2*pi/obj.constituentsV(c) * hours(d_secs));
            end
            Mw = ones(numel(d_time), npars(3));

            for c = 1:numel(obj.constituentsW)
                Mw(:,2*c) = cos(2*pi/obj.constituentsW(c) * hours(d_secs));
                Mw(:,2*c + 1) = sin(2*pi/obj.constituentsW(c) * hours(d_secs));
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
