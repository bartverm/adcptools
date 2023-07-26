classdef TaylorBased < regularization.Regularization
% Regularization requiring a Taylor model
    properties(SetAccess=protected, Dependent)
        min_order (5,:) double {mustBeFinite, mustBeReal, mustBeInteger}
        min_time_order
        min_s_order
        min_n_order
        min_z_order
        min_sig_order
    end
    methods
        function val = get.min_order(obj)
            val = obj.get_min_order();
        end
        function val = get.min_time_order(obj)
            val = obj.min_order(1,:);
        end
        function val = get.min_s_order(obj)
            val = obj.min_order(2,:);
        end
        function val = get.min_n_order(obj)
            val = obj.min_order(3,:);
        end
        function val = get.min_z_order(obj)
            val = obj.min_order(4,:);
        end
        function val = get.min_sig_order(obj)
            val = obj.min_order(5,:);
        end
        function val = find_par(obj,...
                order_in,... % expansion order
                comp_in,... % component
                vars_in,... % expansion variable
                varargin)
            all_ord = [0 1 2];
            all_comps = obj.model.component_names;
            all_vars = obj.model.var_names;

            if nargin < 2 || isempty(order_in)
                order_in = all_ord;
            end
            if nargin < 3 || isempty(comp_in)
                comp_in = all_comps;
            end
            if nargin < 4 || isempty(vars_in)
                vars_in = all_vars;
            end
            % find number of paramter
            val = obj.model.find_par(order_in, comp_in, vars_in,...
                varargin{:});
            val = repmat(val, [obj.mesh.ncells,1]);
        end
    end
    methods(Access = protected)
        function assemble_matrix_private(obj,varargin)
            assert(obj.model_is_taylor,...
                    "Regularization requires a TaylorModel");
            assert(all(obj.model.lumped_orders > obj.min_order, 'all'),...
                ['Taylor model does not meet minimum required ',...
                 'expansion orders. Check the object''s min_order',...
                 'properties.'])
        end
        function tf = model_is_taylor(obj)
            tf = isa(obj.model,"TaylorModel");
        end
    end
    methods(Static, Access = protected, Abstract)
        get_min_order()
    end
end