classdef DataBlock_Traits < matlab.mixin.Heterogeneous & handle
    properties
        name=char.empty;
        header=char.empty;
        postaction=function_handle.empty;
        fields=PD0.DataField_Traits.empty;
    end
    methods
        function field=get_field_by_name(obj,name)
            f_field=strcmp(name,{obj.fields(:).name});
            if any(f_field)
                field=obj.fields(f_field);
            else
                field=PD0.DataField_Traits.empty;
            end
        end
        function set.name(obj,val)
            validateattributes(val, {'char'},{'scalartext'},'set.name','name',2);
            obj.name=val;
        end
        function set.header(obj,val)
            validateattributes(val, {'char'},{'scalartext'},'set.header','header',2);
            assert(regexpi(val,'^[0-9a-f]*$')==1,'set.header:not_hex','header must be a hexadecimal value');
            obj.header=val;
        end
        function set.postaction(obj,val)
            validateattributes(val,{'function_handle'},{'scalar'},'set.postaction','postaction', 2);
            obj.postaction=val;
        end
        function set.fields(obj,val)
            validateattributes(val,{'PD0.DataField_Traits'},{'vector'},'set.fields','fields',2);
            obj.fields=val;
        end
        
    end
end