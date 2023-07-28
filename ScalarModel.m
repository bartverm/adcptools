classdef ScalarModel < DataModel
    properties
        scalar_name(1,:) char = 'eta'
    end
    methods(Access = protected)
        function val = get_component_names(obj)
            val = {obj.scalar_name};
        end
        function val = get_ncomponents(~)
            val = 1;
        end
    end
end 