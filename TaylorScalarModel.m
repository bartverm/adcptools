classdef TaylorScalarModel < TaylorModel & ScalarModel
    methods
        function obj = TaylorScalarModel(varargin)
            obj = obj@TaylorModel(varargin{:});
            obj = obj@ScalarModel(varargin{:});
            order_unass = intersect(obj.unassigned_properties, { ...
                's_order', 'n_order', 'time_order', 'z_order',...
                'sigma_order'});
            for co = 1:numel(order_unass)
                obj.assign_property(order_unass{co}, 0);
            end
        end
        function val = get_model(obj, varargin)
            val = obj.get_model@TaylorModel(varargin{:});
        end
    end
    methods(Access = protected)
        function val = get_names(obj)
            val = obj.get_names@TaylorModel();
        end
        function val = get_ncomponents(~)
            val = 1;
        end
        function val = get_npars(obj)
            val = obj.get_npars@TaylorModel();
        end
    end
end