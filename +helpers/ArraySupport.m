classdef ArraySupport < handle & matlab.mixin.Copyable
% Class to help onstruct object arrays based on array input to constructor
%
%   obj = ArraySupport(...) constructs obj. The size of any non-scalar
%   cell or handle class input determines the size of the output object
%   array. All non scalar inputs of type cell or handle class will need to
%   have the same size. To use this constructor subclass from
%   ArraySupport and call the superclass constructor in your class.
%
%   ArraySupport methods (protected):
%   assign_property - Assign elements of an array to property of obj array
%   run_method - Run method of array of objects
%   
%
%   Example construction:
%   classdef MyClass < helpers.ArraySupport
%       methods
%           function obj = MyClass(varargin)
%               obj = obj@helpers.ArraySupport(varargin{:});
%           end
%  

    methods
        function obj = ArraySupport(varargin)
            if isempty(varargin)
                return
            end
            % figure out which are the non scalar expanded input
            expand_var = cellfun(@(x) (iscell(x) || isa(x,'handle') )...
                && ~isscalar(x), varargin);
            if any(expand_var)
                % get size of expanded variables
                siz_nscal = cellfun(@size,varargin(expand_var), ...
                    'UniformOutput', false);
                % make sure their sizes match
                assert(isscalar(siz_nscal) || isequal(siz_nscal{:}),...
                    'Solver:NonMatchingInputSize',...
                    'Size of all non-scalar input should match')
                siz_obj = num2cell(siz_nscal{1});
                obj(siz_obj{:})=copy(obj);
                for co = 2:numel(obj)-1 % needed for heterogeneous classes
                    obj(co) = copy(obj(1));
                end
            end
        end
    end
    methods(Access = protected, Sealed)
        function assign_property(obj, var_name, var)
% Assign elements of an array to property of object array
%
%   obj.assign_var(var_name, var) assign the variable var to the property 
%   'var_name' of obj. If var is scalar, var is assigned to the propery of
%   each object in the array. If var matches the size of the object array
%   each element in var is assigned to the property of the corresponding
%   object in the array. 
%   If var is a cell, the cell content is assigned and not the cell itself.
%   
%   see also:
%   ArraySupport, run_method
            if isscalar(obj)
                obj.(var_name) = var;
                return
            end
            assert(isscalar(var) || isequal(size(var), size(obj)),...
                'ArrayOnConstruct:WrongInputSize',...
                 ['Input variable should either be scalar or match the',...
                 'size of object array.'])

            if ~iscell(var)
                var=num2cell(var);
            end
            
            [obj.(var_name)] = deal(var{:});
        end

        function varargout = run_method(obj, method_name, varargin)
%   Run a method for all objects in object array obj. 
%
%   [...] = obj.run_method(method_name, ...) runs the method
%   specified by method_name for each object in the object array obj. The
%   input arguments should have the same size as the object array. Each
%   input element is passed to the function and the outputs returned are
%   cells holding the output of the function. The cells have the same size
%   as the object array. If all the outputs of the function are scalars, an
%   array is produced instead of a cell
%
%   See also:
%   ArraySupport, assign_paramater
                if isscalar(obj)
                    feval(method_name, obj, varargin{:});
                    return
                end
                siz_out = size(obj);
                argout = cell(numel(obj),nargout);
                if ~isempty(varargin)
                    siz_in = cellfun(@size, varargin, ...
                        "UniformOutput",false);
                    assert(isequal(size(obj), siz_in{:}),...
                        'helper.ArraySupport:NonMatchingInputSize',...
                        'Size of inputs must match size of object array')
                    varargin = cellfun(@(x) reshape(x,[],1), varargin, ...
                        'UniformOutput',false);
                    not_cell = ~cellfun(@iscell,varargin);
                    varargin(not_cell) = cellfun(@num2cell, ...
                        varargin(not_cell),'UniformOutput',false);
                    varargin = [varargin{:}];
                    
                else 
                    varargin = cell.empty(numel(obj),0);
                end
                for co = 1:numel(obj)
                    [argout{co,:}] = feval(method_name, obj(co), varargin{co,:});
                end

                argout = num2cell(argout,1);
                varargout = cellfun(@(x) obj.reshape_output(x,siz_out), ...
                    argout, 'UniformOutput',false);
        end
    end
    methods
        function ax = generate_tiles(obj)
            for ct = 1:numel(obj)
                ax(ct) = nexttile; %#ok<AGROW> 
            end
        end
    end
    methods(Access = protected, Static)
        function arg = reshape_output(arg,siz_out)
            if all(cellfun(@isscalar,arg))
                arg = [arg{:}];
            end
            arg = reshape(arg,siz_out);
        end
    end
end
