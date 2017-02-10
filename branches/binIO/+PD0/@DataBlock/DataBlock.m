classdef DataBlock < dynamicprops
    properties(SetObservable)
        traits=PD0.DataBlock_Traits.empty(0,0);
    end
    properties(Access=protected)
        dyn_props_meta=meta.DynamicProperty.empty(0,0);
    end
    methods
        function obj=DataBlock(val)
            addlistener(obj,'traits','PostSet',@obj.update_properties);
            if nargin > 0
                obj.traits=val;
            end
        end
        function set.traits(obj,val)
            validateattributes(val,{'PD0.DataBlock_Traits'},{'scalar'},'set.traits','traits',2)
            obj.traits=val;
        end
    end
    methods(Access=protected)
        function update_properties(obj,~,~)
            delete(obj.dyn_props_meta)
            for cnt_field=1:numel(obj.traits.fields)
                obj.dyn_props_meta(cnt_field)=obj.addprop(obj.traits.fields(cnt_field).name);
            end
        end
    end
end