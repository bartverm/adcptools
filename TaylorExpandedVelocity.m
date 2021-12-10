classdef TaylorExpandedVelocity < handle
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