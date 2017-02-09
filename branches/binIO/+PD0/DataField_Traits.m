classdef DataField_Traits < handle
    properties
        name
        size
        type
        offset
    end
    methods
        function set.name(obj,val)
            validateattributes(val, {'char'},{'scalartext'},'set.name','name',2);
            obj.name=val;
        end
        function set.size(obj,val)
            validateattributes(val, {'double'},{'scalar','integer','positive','nonzero'},'set.size','size',2);
            obj.size=val;
        end
        function set.type(obj,val)
            validateattributes(val, {'char'},{'scalartext'},'set.name','name',2);
            assert(any(strcmp(val, {'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', 'uint64', 'int64', 'single', 'double', 'char'})),'set.type:not_valid', 'type must be one of: ''uint8'', ''int8'', ''uint16'', ''int16'', ''uint32'', ''int32'', ''uint64'', ''int64'', ''single'', ''double'', ''char''');
            obj.type=val;
        end
        function set.offset(obj,val)
            validateattributes(val, {'double'},{'scalar','integer','positive','nonzero'},'set.offset','offset',2);
            obj.offset=val;
        end
    end
    
end