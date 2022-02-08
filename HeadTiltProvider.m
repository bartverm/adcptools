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
        function val = heading(obj,adcp)
            if isscalar(obj)
                val = obj.get_heading(adcp);
                return
            end
            fgood = find(obj.has_data,1,'first');
            if isempty(fgood)
                error('No heading availale')
            end
            val = obj(fgood).heading(adcp);
        end
        function val = pitch(obj,adcp)
            if isscalar(obj)
                val = obj.get_pitch(adcp);
                return
            end
            fgood = find(obj.has_data,1,'first');
            if isempty(fgood)
                error('No pitch availale')
            end
            val = obj(fgood).pitch(adcp);
        end
        function val = roll(obj,adcp)
            if isscalar(obj)
                val = obj.get_roll(adcp);
                return
            end
            fgood = find(obj.has_data,1,'first');
            if isempty(fgood)
                error('No roll availale')
            end
            val = obj(fgood).roll(adcp);
        end
    end
    methods(Abstract, Access = protected)
        get_has_data(obj, adcp)
        get_head_tilt_matrix(obj, adcp)
        get_pitch(obj, adcp)
        get_roll(obj, adcp)
        get_heading(obj, adcp)
    end
end