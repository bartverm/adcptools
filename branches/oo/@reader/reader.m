classdef reader < handle
    % Abstract function defining interface of adcp reading functions
    %    Any function that reads adcp data must be subclassed from this
    %    function.
    properties(SetAccess=protected)
        adcp;
        status=reader_status.Invalid;
    end
    methods(Abstract)
        read(obj);
    end
    methods
        function set.adcp(obj,val)
            assert(isa(val,'adcp'));
            obj.adcp=val;
        end
        function set.status(obj,val)
            assert(isa(val,'reader_status'))
            obj.status=val;
        end
    end
end