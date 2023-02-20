classdef TaylorTidalVelocityModel < TidalModel & TaylorModel
    %   Velocity model based on Taylor expansions of tidal amplitudes.
    %
    %   @TaylorTidalVelocityModel
    %   npars - number of parameters to be estimated    
    %
    %   TaylorTidalVelocityModel properties:
    %   @TaylorVelocityModel
    %   s_order - expand to given orders in s coordinate
    %   n_order - expand to given orders in n coordinate
    %   z_order - expand to given orders in z coordinate (not recommended)
    %   sigma_order - expand to given orders in sigma coordinate
    %
    %   @TidalVelocityModel
    %   constituentsU - names of constituents of velocity in x/s direction
    %   constituentsV - names of constituents of velocity in y/n direction
    %   constituentsW - names of constituents of velocity in z direction
    %
    %   TaylorTidalVelocityModel methods:
    %
    %   [Mu, Mv, Mw]=get_model(obj, time, d_s, d_n, d_z, d_sigma)
    %   procuces the model matrices consisting of Kronecker
    %   products of the Taylor and Tidal velocity models,
    %   respectively.
    %
    %   [p, pu, pv, pw] = get_parameter_names(obj)
    %   produces a cell of strings of length npars containing the names of
    %   the parameters that are estimated. p = [pu pv pw].

    %  see also: VelocityModel, TaylorVelocityModel, TidalVelocityModel

    properties
%         % TaylorExpandedVelocity/s-order order of expansion in s-direction
%         %
%         %   1x3 row vector holding non-negative integers that specify up
%         %   to which order the velocity should be expanded in the
%         %   s-direction for each of the three cartesian velocity components
%         %   [2 1 0] means expand x-velocity up to second order in
%         %   s-direction, y-velocity up to first order and z-velocity up to
%         %   0th order. For the x-velocity this results in three model
%         %   parameters (mean x-vel, dx-vel/ds, d^2 x_vel/ds^2), for the
%         %   y-velocity in two parameters (mean y-vel, dy-vel/dx) and for
%         %   z-velocity in 1 parameter (mean z-vel).
%         %
%         %   Default is: [0, 0, 0]
%         %
%         %   see also: TaylorExpandedVelocity, n_order, z_order,
%         %   sigma_order, time_order
%         s_order (1,3) double {...
%             mustBeFinite, mustBeInteger, mustBeNonnegative} = [0 0 0];
% 
%         % TaylorExpandedVelocity/n-order order of expansion in n-direction
%         %
%         %   1x3 row vector holding non-negative integers that specify up
%         %   to which order the velocity should be expanded in the
%         %   n-direction for each of the three cartesian velocity
%         %   components.
%         %
%         %   Default is: [0, 0, 0]
%         %
%         %   See for an example: s_order
%         %
%         %   see also: TaylorExpandedVelocity, s_order, z_order,
%         %   sigma_order, time_order
%         n_order (1,3) double {...
%             mustBeFinite, mustBeInteger, mustBeNonnegative} = [0 0 0];
% 
%         % TaylorExpandedVelocity/z-order order of expansion in z-direction
%         %
%         %   1x3 row vector holding non-negative integers that specify up
%         %   to which order the velocity should be expanded in the
%         %   z-direction for each of the three cartesian velocity
%         %   components.
%         %
%         %   Default is: [0, 0, 0]
%         %
%         %   See for an example: s_order
%         %
%         %   see also: TaylorExpandedVelocity, s_order, n_order,
%         %   sigma_order, time_order
%         z_order (1,3) double {...
%             mustBeFinite, mustBeInteger, mustBeNonnegative} = [0 0 0];
% 
%         % TaylorExpandedVelocity/sigma-order order of expansion in sigma direction
%         %
%         %   1x3 row vector holding non-negative integers that specify up
%         %   to which order the velocity should be expanded in the
%         %   sigma-direction for each of the three cartesian velocity
%         %   components.
%         %
%         %   Default is: [0, 0, 0]
%         %
%         %   See for an example: s_order
%         %
%         %   see also: TaylorExpandedVelocity, s_order, z_order,
%         %   n_order, time_order
%         sigma_order (1,3) double {...
%             mustBeFinite, mustBeInteger, mustBeNonnegative} = [0 0 0];
% 
%         % TaylorExpandedVelocity/time-order order of expansion in time
%         %
%         %   1x3 row vector holding non-negative integers that specify up
%         %   to which order the velocity should be expanded in
%         %   time for each of the three cartesian velocity
%         %   components.
%         %
%         %   Default is: [0, 0, 0]
%         %
%         %   See for an example: s_order
%         %
%         %   see also: TaylorExpandedVelocity, s_order, z_order,
%         %   n_order, n_order
%         %   time_order (1,3) double {...
%         % mustBeFinite, mustBeInteger, mustBeNonnegative} = [0 0 0];
% 
%         % TidalVelocityModel/constituentsU constituents for x-component of velocity
%         %
%         %   1xN row vector defining tidal constituents to be included in the model
%         %   for the x-component of the velocity. Every value given indicates the
%         %   period of the constituent to be included. For every given value, two
%         %   model parameters are fitted, which are the coefficients for the sin and
%         %   cos function. From those amplitude and phases are computed. A residual
%         %   is always included.
%         %
%         %   see also: TidalVelocityModel, constituentsV, constituentsW,
%         %               get_tidal_pars
%         constituentsU(1,:) cell = {};
% 
%         % TidalVelocityModel/constituentsV constituents for y-component of velocity
%         %
%         %   1xN row vector defining tidal constituents to be included in the model
%         %   for the y-component of the velocity. Every value given indicates the
%         %   period of the constituent to be included. For every given value, two
%         %   model parameters are fitted, which are the coefficients for the sin and
%         %   cos function. From those amplitude and phases are computed. A residual
%         %   is always included.
%         %
%         %   see also: TidalVelocityModel, constituentsU, constituentsW,
%         %               get_tidal_pars
%         constituentsV(1,:) cell = {};
% 
%         % TidalVelocityModel/constituentsW constituents for z-component of velocity
%         %
%         %   1xN row vector defining tidal constituents to be included in the model
%         %   for the z-component of the velocity. Every value given indicates the
%         %   period of the constituent to be included. For every given value, two
%         %   model parameters are fitted, which are the coefficients for the sin and
%         %   cos function. From those amplitude and phases are computed. A residual
%         %   is always included.
%         %
%         %   see also: TidalVelocityModel, constituentsU, constituentsV,
%         %               get_tidal_pars
%         constituentsW(1,:) cell = {};


%         TaylorTidalVelocityModel/names names of the different parameter
%         names. The are coded as

          names(1,:) cell = {};
    end

    methods
        function [Mu, Mv, Mw]=get_model(obj, time, d_s, d_n, d_z, d_sigma)
            %   [Mu, Mv, Mw]=get_model(obj, time, d_s, d_n, d_z, d_sigma)
            %   procuces the model matrices consisting of Kronecker
            %   products of the Taylor and Tidal velocity models,
            %   respectively. 
            [Tu, Tv, Tw] = get_model@TidalVelocityModel(obj,time);
            [Su, Sv, Sw] = get_model@TaylorVelocityModel(obj, time, d_s, d_n, d_z, d_sigma);
            ns = numel(time);
            Mu = zeros([ns,obj.npars(1)]);
            Mv = zeros([ns,obj.npars(2)]);
            Mw = zeros([ns,obj.npars(3)]);

            for i = 1:ns
                Mu(i,:) = kron(Su(i,:), Tu(i,:));
                Mv(i,:) = kron(Sv(i,:), Tv(i,:));
                Mw(i,:) = kron(Sw(i,:), Tw(i,:));
            end

        end

        function obj = get_parameter_names(obj)
            pu = tidal_columns(obj, 'u'); npu = length(pu);
            pv = tidal_columns(obj, 'v'); npv = length(pv);
            pw = tidal_columns(obj, 'w'); npw = length(pw);
            p = cell([1, sum(obj.npars)]);
            tid_names = {pu, pv, pw};
            np = {npu, npv, npw};
            varnames = {'u', 'v', 'w'};

            idx = 1;
            for vari = 1:3
                for i=1:np{vari}
                    p{idx} = obj.tidal_taylor_name(varnames{vari}, 0, 0, tid_names{vari}{i});
                    idx = idx + 1;
                end
                for ord = 1:obj.s_order(vari)
                    for i = 1:np{vari}
                        p{idx} = obj.tidal_taylor_name(varnames{vari}, 'x', ord, tid_names{vari}{i});
                        idx = idx + 1;
                    end
                end
                for ord = 1:obj.n_order(vari)
                    for i = 1:np{vari} %
                        p{idx} = obj.tidal_taylor_name(varnames{vari}, 'y', ord, tid_names{vari}{i});
                        idx = idx + 1;
                    end
                end

                for ord = 1:obj.z_order(vari)
                    for i = 1:np{vari}
                        p{idx} = obj.tidal_taylor_name(varnames{vari}, 'z', ord, tid_names{vari}{i});
                        idx = idx + 1;
                    end
                end

                for ord = 1:obj.sigma_order(vari)
                    for i = 1:np{vari}
                        p{idx} = obj.tidal_taylor_name(varnames{vari}, 'sig', ord, tid_names{vari}{i});
                        idx = idx + 1;
                    end
                end
            end
%             pu = p(1:obj.npars(1));
%             pv = p(obj.npars(1)+1:obj.npars(1) + obj.npars(2));
%             pw = p(obj.npars(1) + obj.npars(2) + 1:end);
            obj.names = p;
        end

        function p = tidal_columns(obj, var)
            switch var
                case 'u'
                    constituents = obj.constituentsU;
                case 'v'
                    constituents = obj.constituentsV;
                case 'w'
                    constituents = obj.constituentsW;
            end
            p = cell([1, 2*numel(constituents)+1]);
            p{1} = 'M0A';
            for j = 1:numel(constituents)
                p{2*j} = sprintf('%s%s', constituents{j}, 'A');
                p{2*j+1} = sprintf('%s%s', constituents{j}, 'B');
            end
        end

        function par_name = tidal_taylor_name(obj, var, dim, ord, tid_name)
            tay_name = obj.taylor_name(var, dim, ord);
            par_name = sprintf('%s: %s', tay_name, tid_name);
        end

        function par_name = taylor_name(obj, var, dim, ord)
            if ord == 0
                par_name = sprintf('%s', [var, '0']);
            else
                par_name = sprintf('d^%u%s/d%s^%u', ord, var, dim, ord);
            end
        end


    end

    methods(Access=protected)
        function val = get_npars(obj)
            val = ones(1,3) +...
                obj.s_order +...
                obj.n_order +...
                obj.z_order +...
                obj.sigma_order;

            val = val.*[2*numel(obj.constituentsU) + 1, ...
                2*numel(obj.constituentsV) + 1, ...
                2*numel(obj.constituentsW) + 1];
        end
        function val = get_ncomponents(obj)

        end
    end
end


