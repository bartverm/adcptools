classdef VelocityModel < DataModel &...
        matlab.mixin.Heterogeneous
    methods
        function M = get_model(obj, d_time, ~, ~, ~, ~)
            M = ones(numel(d_time), 1, obj.ncomponents);
        end
    end
    methods(Access=protected)
        function val = get_ncomponents(~)
            val = 3;
        end
        function val = get_component_names(~)
            val = {'u', 'v', 'w'};
        end
        function val=get_npars(~)
            val = [1 1 1];
        end
        function val=get_names(~)
            val = {{'u'}, {'v'}, {'w'}};
        end
    end
end