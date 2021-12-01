classdef HeadingProvider < matlab.mixin.Heterogeneous & handle
    methods (Sealed)
        function val=has_data(obj,adcp)
            if isscalar(obj)
                val=obj.get_has_data(adcp);
            else
                val=false(size(obj));
                for cel=1:numel(obj)
                    val(cel)=obj(cel).get_has_data(adcp);
                end
            end
        end
        function val=heading(obj,adcp)
            if isscalar(obj)
                val=obj.get_heading(adcp);
            else
                idx=find(obj.has_data(adcp),1,"first");
                val=obj(idx).get_heading(adcp);
            end
        end
    end
    methods (Abstract, Access=protected)
        val=get_heading(obj,adcp)
        val=get_has_data(obj,adcp)
    end
end