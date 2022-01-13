classdef HeadTiltProvider < handle & matlab.mixin.Heterogeneous
    methods(Sealed)
        function val = has_data(obj, adcp)
            if isscalar(obj)
                val = obj.get_has_data(adcp);
                return
            end
            val = false(size(obj));
            for co = 1:numel(obj)
                val(co) = obj(co).has_data(adcp);
            end
        end
        function val = head_tilt_matrix(obj, adcp)
            if isscalar(obj)
                val = obj.get_head_tilt_matrix(adcp);
                return
            end
            fgood = find(obj.has_data,1,'first');
            if isempty(fgood)
                error('No heading and tilt matrix availale')
            end
            val = obj(fgood).head_tilt_matrix(adcp);
        end
    end
    methods(Abstract, Access = protected)
        get_has_data(obj, adcp)
        get_head_tilt_matrix(obj, adcp)
    end
end