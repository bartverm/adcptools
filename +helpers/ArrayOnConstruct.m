classdef ArrayOnConstruct < handle
    methods
        function obj = ArrayOnConstruct(varargin)
            if isempty(varargin)
                return
            end
            expand_var = cellfun(@(x) (iscell(x) || isa(x,'handle') )...
                && ~isscalar(x), varargin);
            if any(expand_var)
                siz_nscal = cellfun(@size,varargin(expand_var),'UniformOutput', false);
                assert(isscalar(siz_nscal) || isequal(siz_nscal{:}),...
                    'Solver:NonMatchingInputSize',...
                    'Size of all non-scalar input should match')
                siz_obj = num2cell(siz_nscal{1});
                obj(siz_obj{:})=obj;
            end

        end
    end
    methods(Access = protected)
        function assign_var(obj, var_name, var)
            if ~iscell(var) && isequal(size(obj),size(var))
                var = num2cell(var);
            elseif ~iscell(var)
                var = {var};
            end
            assert(isscalar(var) || isequal(size(obj), size(var)), ...
                'Solver:NonMatchingInputSize',...
                'Size of all non-scalar input should match with object array size')
            [obj.(var_name)] = deal(var{:});
        end
    end
end