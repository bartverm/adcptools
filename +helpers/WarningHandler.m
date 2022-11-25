classdef WarningHandler < handle
% helper class to handle warnings
%
%   WarningHandler properties:
%   disabled_warnings - disabled warnings and their original state
%
%   WarningHandler methods:
%   warn_and_disable - issue warning and disable afterwards
%
%   Upon deletion of the object, disabled warning will be restored to their
%   original state.
    properties(GetAccess = public, SetAccess = private)
        disabled_warnings(:,1) struct = struct('identifier', {},...
            'state', {})
    end
    methods
        function delete(obj)
            warning(obj.disabled_warnings)
        end
        function warn_and_disable(obj, warn_id, warn_message)
            if any(strcmp(warn_id,...
                    {obj.disabled_warnings.identifier}))
                return
            end
            warning(warn_id, [warn_message,...
                '\nWarning is disabled until originating',...
                'object is deleted.'])
            obj.disabled_warnings = [obj.disabled_warnings;...
                warning('off', warn_id)];
        end
    end
end