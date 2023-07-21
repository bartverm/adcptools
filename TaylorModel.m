classdef TaylorModel < DataModel
    % Velocity model based on Taylor expansions.
    %
    %   TaylorModel properties:
    %   time_order - expand to given orders in time
    %   s_order - expand to given orders in s coordinate
    %   n_order - expand to given orders in n coordinate
    %   z_order - expand to given orders in z coordinate
    %   sigma_order - expand to given orders in sigma coordinate
    %
    %   TaylorModel methods:
    %
    %  see also: VelocityModel, DataModel
    properties(Constant)
        var_names = {'t', 's', 'n', 'z', 'sig'};
    end

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
        npars_not_expanded
        lumped_orders
    end

    methods
        function val = get.npars_not_expanded(obj)
            val = obj.npars ./ obj.get_npars_tay;
        end
        function val = get.lumped_orders(obj)
            val = obj.lump_orders;
        end
        function set.time_order(obj,val)
            obj.check_order(val);
            obj.time_order = val;
        end
        function set.s_order(obj,val)
            obj.check_order(val);
            obj.s_order = val;
        end
        function set.n_order(obj,val)
            obj.check_order(val);
            obj.n_order = val;
        end
        function set.z_order(obj,val)
            obj.check_order(val);
            obj.z_order = val;
        end
        function set.sigma_order(obj,val)
            obj.check_order(val);
            obj.sigma_order = val;
        end
        function val = find_par(obj,...
                order_in,... % expansion order
                comp_in,... % component
                vars_in)    % expansion variable
            % select given parameters
            %
            % idx = obj.find_par(order) find index of terms with the given
            %   order (0, 1 or 2). If it is empty, returns terms of all
            %   orders
            %
            % idx = obj.find_par(order, comp) find index of terms of given
            %   order and component (one of obj.component_names). If empty
            %   returs all 
            %
            % idx = obj.find_par(order, comp, var) find index of terms of
            %   given order, component, and expanded to the given variable.
            %   For second order, all double and cross-derivatives with the
            %   given variable are returned
            all_comps = obj.get_component_names;
            all_vars = obj.var_names;
            if nargin < 2
                val = true(sum(obj.npars),1);
                return
            end
            if isempty(order_in)
                order_in = [0 1 2];
            end
            if nargin < 3 || isempty(comp_in)
                comp_in = all_comps;
            end
            if nargin < 4 || isempty(vars_in)
                vars_in = all_vars;
            end
            assert(iscellstr(comp_in) || ischar(comp_in),...
                'Components should be given as char or cell of chars')
            if ischar(comp_in)
                comp_in = {comp_in};
            end
            assert(all(ismember(comp_in, all_comps)),...
                ['Components should be one or more of ''',...
                strjoin(all_comps),...
                ''''])
            assert(iscellstr(vars_in) || ischar(vars_in),...
                'Expansion variable should be a char or a cell of chars')
            if ischar(vars_in)
                vars_in = {vars_in};
            end
            assert(all(ismember(vars_in, all_vars)),...
                ['Expansion variables should be one or more of ''',...
                strjoin(all_vars),...
                ''''])
            assert(isnumeric(order_in),'Order should be numeric')
            assert(isreal(order_in), 'Order should be a real number')
            assert(all(isfinite(order_in)), 'Order should be finite')
            assert(all(order_in >= 0, 'all'), 'Order should not be negative')
            assert(all(order_in < 3, 'all'), 'Order cannot be larger than two');
            assert(all(floor(order_in)==order_in,'all'),...
                'Order should contain only integer numbers')
            lo = obj.lump_orders();
            np = sum(obj.get_npars_tay);
            nc = obj.ncomponents;
            vnames = obj.var_names;
            val = false(np,1);
           
            pos = 1;
            for comp = 1:nc
                comp_name = all_comps(comp);

                % order 0
                ord = 0;
                if ismember(ord, order_in) &&...
                        ismember(comp_name, comp_in)
                    val(pos) = true;
                end
                pos = pos + 1;
                
                % order 1
                ord = 1;
                exp_vars = find(lo(:,comp)>0);
                for cvar = 1:numel(exp_vars)
                    cur_var = exp_vars(cvar);
                    var_name = vnames(cur_var);
                    if ismember(ord, order_in) &&...
                            ismember(comp_name, comp_in) &&...
                            ismember(var_name, vars_in)
                        val(pos) = true;
                    end
                    pos = pos + 1;
                end

                % order 2
                ord = 2;
                exp_vars = find(lo(:,comp)>1);
                for cvar1 = 1:numel(exp_vars)
                    cur_var1 = exp_vars(cvar1);
                    var_name1 = vnames(cur_var1);
                    if ismember(ord, order_in) &&...
                            ismember(comp_name, comp_in) &&...
                            ismember(var_name1, vars_in)
                        val(pos) = true;
                    end
                    pos = pos + 1;
                    for cvar2 = cvar1+1:numel(exp_vars)
                        cur_var2 = exp_vars(cvar2);
                        var_name2 = vnames(cur_var2);
                        if ismember(ord, order_in) &&...
                                ismember(comp_name, comp_in) &&...
                                ismember(var_name1, vars_in) &&...
                                ismember(var_name2, vars_in)
                            val(pos) = true;
                        end
                        pos = pos + 1;
                    end                     
                end                
            end
            
        end
        function val = get.npars_per_order(obj)
            lo = obj.lump_orders();
            nc = obj.ncomponents;
            no1 = sum(lo>0,1); % number of order 1 parameters
            no2 = sum(lo>1,1); % number of order 2 parameter
            val = [ones(1,nc);... % zero order always included
                no1;...        % first order terms
                (no2.^2+no2)/2]; % second order independent terms: 1+2+...+N = (N^2 + N)/2
        end
        function M = get_model(obj, time, d_s, d_n, d_z, d_sigma)
            np_per_ord = obj.npars_per_order;
            max_np = max(obj.get_npars_tay);
            nc = obj.ncomponents;
            lo = obj.lump_orders;
            % Convert dt to seconds: numerical values
            d_t = seconds(time - mean(time));

            M = nan(numel(time), max_np, nc);
            dev_vec = [d_t, d_s, d_n, d_z, d_sigma]; % deviatoric vector (x - a) for first order terms
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
    end
    methods(Static)
        function par_name = tayl_name(comp, coord, ord)
            if ord == 0
                par_name = sprintf('%s', [comp, '0']);
            elseif ord == 1
                par_name = sprintf('d%s/d%s', comp, coord);
            elseif ord == 2
                if iscell(coord)
                    par_name = sprintf('d^2%s/d%sd%s', comp, coord{1},...
                        coord{2});
                else
                    par_name = sprintf('d^2%s/d%s^2',comp, coord);
                end
            end
        end

    end
    methods (Access = protected)
        function check_order(obj,val)
            assert(numel(val) == obj.ncomponents, 'The number of elements should match the number of components')
        end
        function names = get_names(obj)
            % Forms a cell array of dimensions 1xobj.ncomponents
            % Elements of the cell array are cell arrays containing the
            % Taylor expansion names per component
            coord_names = obj.var_names;
            names = cell([obj.ncomponents, 1]);
            lo = obj.lump_orders;
            ncomp = obj.ncomponents;
            npars = obj.get_npars_tay; % calling protected method for subclasses
            for comp = 1:ncomp
                comp_name = obj.component_names{comp};
                names{comp} = cell(1,npars(comp));

                % order zero
                names{comp}{1} = [obj.component_names{comp}, '0']; % Spatially constant part

                % order 1
                idx = 2;
                vars_ord1 = find(lo(:,comp)>0);
                for cv = 1:numel(vars_ord1)
                    cur_var = vars_ord1(cv);
                    var_name = coord_names{cur_var};
                    names{comp}{idx} = obj.tayl_name(comp_name,...
                        var_name, 1);
                    idx = idx + 1;
                end

                % order 2
                vars_ord2 = find(lo(:,comp)>1);
                for cv1 = 1:numel(vars_ord2)
                    cur_var1 = vars_ord2(cv1);
                    var_name1 = coord_names{cur_var1};
                    names{comp}{idx} = obj.tayl_name(comp_name,...
                        var_name1, 2); % diagonal terms
                    idx = idx + 1;
                    for cv2 = cv1+1:numel(vars_ord2)
                        cur_var2 = vars_ord2(cv2);
                        var_name2 = coord_names{cur_var2};
                        names{comp}{idx} = obj.tayl_name(comp_name,...
                            {var_name1, var_name2}, 2); % cross derivatives
                        idx = idx + 1;
                    end
                end
            end
        end
    
        function val = get_npars(obj)
            val = obj.get_npars_tay;
        end

        function val = get_npars_tay(obj)
            val = obj.npars_per_order();
            val = sum(val,1);
        end    

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

    end


end
