classdef ScalarModel < DataModel &...
        matlab.mixin.Heterogeneous
    properties
        scalar_name(1,:) char = 'eta'
    end
    methods
        function M = get_model(~, d_time, ~, ~, ~, ~)
            M = ones(numel(d_time), 1, 1);
        end
    end
    methods(Access = protected)
        function val = get_component_names(obj)
            val = {obj.scalar_name};
        end
        function val = get_ncomponents(~)
            val = 1;
        end
        function val=get_npars(~)
            val = 1;
        end
        function val=get_names(obj)
            val = {{obj.scalar_name}};
        end
    end
end 