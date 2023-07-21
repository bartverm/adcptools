classdef TaylorBasedRegularization < Regularization
% Regularization requiring a Taylor model
    methods
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
            if ~obj.model_is_taylor
                error("TaylorBasedRegularizon:NoTaylorModel",...
                    "Regularization requires a TaylorModel");
            end
        end
        function tf = model_is_taylor(obj)
            tf = isa(obj.model,"TaylorModel");
        end
    end
end