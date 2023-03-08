classdef TaylorModel < DataModel
    % Velocity model based on Taylor expansions.
    %
    %   TaylorModel properties:
    %   s_order - expand to given orders in s coordinate
    %   n_order - expand to given orders in n coordinate
    %   z_order - expand to given orders in z coordinate
    %   sigma_order - expand to given orders in sigma coordinate
    %   time_order - expand to given orders in time
    %
    %   TaylorModel methods:
    %   get_velocity - return velocity based on model parameters
    %
    %  see also: DataModel

    properties
        % TaylorModel/time-order order of expansion in time
        %
        %   1xNCOMP row vector holding non-negative integers that specify up
        %   to which order the velocity should be expanded with respect to
        %   time for each of the NCOMP components
        %   components. Maximum order is 2
        %
        %   Default is: [0, 0, 0]
        %
        %   See for an example: s_order
        %
        %   see also: TaylorModel, s_order, z_order,
        %   n_order, n_order
        time_order (1,:) double {...
            mustBeFinite,...
            mustBeInteger,...
            mustBeNonnegative,...
            mustBeLessThanOrEqual(time_order,2)...
            } = [0 0 0];

        % TaylorModel/s-order order of expansion in s direction
        %
        %   1xNCOMP row vector holding non-negative integers that specify up
        %   to which order the velocity should be expanded with respect to
        %   s (accros cross section) for each of the NCOMP components
        %   components. Maximum order is 2
        %
        %   Example:
        %   If we would specify for s-order [2 1 0] and, n-order [2 2 0]
        %   and for all the others [0 0 0] it means we want to expand the
        %   first component to the second order in the s variable and to
        %   the second order in the n-direction. We get the following
        %   parameters for the first component (say u):
        %   u0 du/ds du/dn d^2u/ds^2 d^2u/dnds d^2u/dn^2 (6 parameters)
        %   For the second component we expand second order to s and only
        %   first order to n. We get the following parameters for the
        %   second component (say v):
        %   v0 dv/ds dv/dn d^2v/ds^2 (4 parameters)
        %   Note the order in which parameters will be estimated!
        %
        %   Default is: [0, 0, 0]
        %
        %
        %   see also: TaylorModel, s_order, z_order,
        %   n_order, n_order
        s_order (1,:) double {...
            mustBeFinite,...
            mustBeInteger,...
            mustBeNonnegative,...
            mustBeLessThanOrEqual(s_order,2)...
            } = [0 0 0];

        % TaylorModel/n-order order of expansion in n direction
        %
        %   1xNCOMP row vector holding non-negative integers that specify up
        %   to which order the velocity should be expanded with respect to
        %   n (along cross section) for each of the NCOMP components
        %   components. Maximum order is 2
        %
        %   Default is: [0, 0, 0]
        %
        %   See for an example: s_order
        %
        %   see also: TaylorModel, s_order, z_order,
        %   n_order, n_order
        n_order (1,:) double {...
            mustBeFinite,...
            mustBeInteger,...
            mustBeNonnegative,...
            mustBeLessThanOrEqual(n_order,2)...
            } = [0 0 0];

        % TaylorModel/z-order order of expansion in z direction
        %
        %   1xNCOMP row vector holding non-negative integers that specify up
        %   to which order the velocity should be expanded with respect to
        %   z (vertical) for each of the NCOMP components
        %   components. Maximum order is 2
        %
        %   Default is: [0, 0, 0]
        %
        %   See for an example: s_order
        %
        %   see also: TaylorModel, s_order, z_order,
        %   n_order, n_order
        z_order (1,:) double {...
            mustBeFinite,...
            mustBeInteger,...
            mustBeNonnegative,...
            mustBeLessThanOrEqual(z_order,2)...
            } = [0 0 0];

        % TaylorModel/sigma-order order of expansion in sigma
        %
        %   1xNCOMP row vector holding non-negative integers that specify up
        %   to which order the velocity should be expanded with respect to
        %   sigma (normalized vertical) for each of the NCOMP components
        %   components. Maximum order is 2
        %
        %   Default is: [0, 0, 0]
        %
        %   See for an example: s_order
        %
        %   see also: TaylorModel, s_order, z_order,
        %   n_order, n_order
        sigma_order (1,:) double {...
            mustBeFinite,...
            mustBeInteger,...
            mustBeNonnegative,...
            mustBeLessThanOrEqual(sigma_order,2)...
            } = [0 0 0];
    end
    properties(Dependent, SetAccess=protected, GetAccess = public)
        npars_per_order
        names
    end
    methods(Access = protected)
        function check_components(obj)
            assert(isequal(...
                numel(obj.time_order),...
                numel(obj.s_order),...
                numel(obj.n_order),...
                numel(obj.z_order),...
                numel(obj.sigma_order)...
                ),...
                'TaylorModel:NonMatchingComponents',...
                'Number of components in each expansion variable should be equal');
        end
        function val = get_ncomponents(obj)
            obj.check_components();
            val = numel(obj.n_order);
        end
        function val = lump_orders(obj)
            obj.check_components();
            val = [obj.time_order;...
                obj.s_order;...
                obj.n_order;...
                obj.z_order;...
                obj.sigma_order];
        end

        function val = get_npars(obj)
            val = obj.npars_per_order();
            val = sum(val,1);
        end
    end
    methods
        function val = get.npars_per_order(obj)
            lo = obj.lump_orders();
            nc = obj.ncomponents;
            no1 = sum(lo>0,1); % number of order 1 parameters
            no2 = sum(lo>1,1); % number of order 2 parameter
            val = [ones(1,nc);... % zero order always included
                no1;...        % first order terms
                (no2.^2+no2)/2]; % second order independent terms: 1+2+...+N = (N^2 + N)/2
        end
        function M=get_model(obj, time, d_s, d_n, d_z, d_sigma)
            np_per_ord = obj.npars_per_order;
            max_np = max(obj.npars);
            nc = obj.ncomponents;
            lo = obj.lump_orders;
            M = nan(numel(time), max_np, nc);
            dev_vec = [time, d_s, d_n, d_z, d_sigma]; % deviatoric vector (x - a) for first order terms
            dev_vec_sq = helpers.matmult(dev_vec,...
                permute(dev_vec,[1 3 2]),2,3); %(x - a)^T(x-a) for second order terms

            pidx = cumsum(np_per_ord,1); % parameter indices in M
            for c_comp = 1:obj.ncomponents
                M(:,1,c_comp) = 1; % zero-th order

                M(:,2:pidx(2,c_comp),c_comp)=...
                    dev_vec(:,lo(:,c_comp)>0); % first order terms

                idx2 = pidx(2,c_comp)+1; % start index for order 2 terms
                n_sec = sum(lo(:,c_comp)>1,1);
                coeff = (2-eye(n_sec))/2; % coefficients (1/2 for diagonal elements, others 1)
                psec = find(lo(:,c_comp)>1); % index of paramaters to expand to second order
                for co = 1:n_sec
                    for ci = co:n_sec
                        M(:,idx2,c_comp) = coeff(co,ci)*...
                            dev_vec_sq(:,psec(co),psec(ci));
                        idx2 = idx2 + 1;
                    end
                end
            end
        end

        function names = get.names(obj)
            % Only supports diagonal Hessian matrices

            % Forms a cell array of dimensions 1xobj.ncomponents
            % Elements of the cell array are cell arrays containing the
            % Taylor expansion names per component
            orders = {obj.time_order, obj.s_order, obj.n_order, obj.z_order, obj.sigma_order};
            coord_names = {'t', 'x', 'y', 'z', 'sig'};
            names = cell([obj.ncomponents, 1]);
            for comp = 1:obj.ncomponents
                names{comp}{1} = [obj.components{comp}, '0']; % Spatially constant part
                idx = 2;
                for coord = 1:length(orders)
                    cur_order = orders{coord}(comp);
                    for ord = 1:cur_order
                        names{comp}{idx} = obj.tayl_name(obj.components{comp}, coord_names{coord}, ord);
                        idx = idx + 1;
                    end
                end
            end
        end

        function par_name = tayl_name(obj, comp, coord, ord)
            if ord == 0
                par_name = sprintf('%s', [comp, '0']);
            else
                par_name = sprintf('d^%u%s/d%s^%u', ord, comp, coord, ord);
            end
        end


    end
end
