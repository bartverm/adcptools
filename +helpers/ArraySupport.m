classdef ArraySupport < handle & matlab.mixin.Copyable
% Class to help onstruct object arrays based on array input to constructor
%
%   obj = ArraySupport(...) constructs obj. The size of any non-scalar
%   cell or handle class input determines the size of the output object
%   array. All non scalar inputs of type cell or handle class will need to
%   have the same size. To use this constructor subclass from
%   ArraySupport and call the superclass constructor in your class.
%
%   obj = ArraySupport('NoExpand',...) any input given after this argument
%   is not used to determine size of constructed array
%
% 
%   ArraySupport methods:
%   assign_property - Assign elements of an array to property of obj array
%   cat_property - Concatenate and return properties from objects
%   run_method - Run method, spreading inputs over obect calls
%   run_method_single - Run method, with same input on object calls
%   
%
%   Example construction:
%   classdef MyClass < helpers.ArraySupport
%       methods
%           function obj = MyClass(varargin)
%               obj = obj@helpers.ArraySupport(varargin{:});
%           end
%  
    properties(Access = protected)
        unprocessed_construction_inputs (1,:) cell = {};
    end
    methods
        function obj = ArraySupport(varargin)
            if isempty(varargin)
                return
            end
            f_expand = find(strcmp('NoExpand',varargin));
            nin = numel(varargin);
            obj(1).unprocessed_construction_inputs = varargin;
            if ~isempty(f_expand) 
                obj(1).unprocessed_construction_inputs(f_expand(1)) = [];
                varargin(f_expand(1):nin) = [];
            end
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
                obj(siz_obj{:}) = copy(obj);
                for co = 2:numel(obj)-1 % needed for heterogeneous classes
                    obj(co) = copy(obj(1));
                end
            end
        end
    end
    methods(Access = public, Sealed)
        function assign_property(obj, var_name, var, varargin)
% Assign elements of an array to property of object array
%
%   obj.assign_var(var_name, var) assign the variable var to the property 
%   'var_name' of obj. If var is scalar, var is assigned to the propery of
%   each object in the array. If var matches the size of the object array
%   each element in var is assigned to the property of the corresponding
%   object in the array. 
%   If var is a cell, the cell content is assigned and not the cell itself.
%   
%   obj.assign_var(...,'NoExpand') the variable is assigned as is to each
%   of the array elements
%
%   'NoExpand'
%   see also:
%   ArraySupport, run_method

            assert(isscalar(obj) || isscalar(var) || isequal(size(var), size(obj)),...
                'ArrayOnConstruct:WrongInputSize',...
                 ['Input variable should either be scalar or match the',...
                 'size of object array.'])

            if isscalar(obj)
                if iscell(var) && isscalar(var)
                    var = var{1};
                end
                obj.(var_name) = var;
                return
            end
            
            if nargin > 3 && strcmp(varargin{1}, 'NoExpand')
                [obj.(var_name)] = deal(var);
            else   
                if ~iscell(var)
                    var=num2cell(var);
                end
                [obj.(var_name)] = deal(var{:});
            end

        end

        function val = cat_property(obj, var_name, varargin)
% Concatenate and return property of objects in array
%
%   val = obj.cat_property(var_name) concatenates the property with
%   name 'var_name' along the 2nd dimension. If sizes along dimension dim
%   are not consistent, output variables are grown along dimension dim to
%   match the variable with the largest size along dim.
%
%   val = obj.cat_property(...,'PropName', PropValue) allows to specify
%   PropertyName, PropertyValue pairs as follows:
%       'FillVal'
%       Specify the value to grow the arrays with. For handle arrays, the 
%       object is grown with unique new objects. If 'FillVal' is specified,
%       it will be filled with copies of 'FillVal'
%
%       'Dimension'
%       Modify the dimension to concatenate along.
%
%   val = obj.cat_property(...) if var_name is a method, any other input
%       passed to function is passed on the method.
%
            % Handle FillVal argument to be passed to grow_array function
            ga_args={};
            fval_pos = find(strcmp(varargin,'FillVal'));
            if ~isempty(fval_pos) && numel(varargin) > fval_pos(1)
                ga_args = varargin(fval_pos(1)+[0 1]);
                varargin(fval_pos(1)+[0 1]) = [];
            end
            
            % Handle Dimension
            dim = 2;
            fval_dim = find(strcmp(varargin,'Dimension'));
            if ~isempty(fval_dim) && numel(varargin) > fval_dim(1)
                dim = varargin(fval_dim(1)+1);
            end


            % scalar case
            if isscalar(obj)
                if isprop(obj,var_name)
                    val = obj.(var_name);
                else
                    val = obj.(var_name)(varargin{:});
                end
                return
            end

            % differentiate call for property or method
            if isprop(obj, var_name)
                val = {obj.(var_name)};
            else
                val = cell(size(obj));
                [val{:}] = obj.(var_name)(varargin{:});
            end

            % grow and concatenate ouput
            [val{:}] = helpers.grow_array(val{:},'SkipDims',dim, ga_args{:});
            val = cat(dim, val{:});
        end

        function out = run_method_single(obj, narg, method_name, varargin)
%   run method for each object in array and return in cell with outputs
%
%   Only works for single output methods! 
% 
%  returns a cell holding each of the outputs of 'method_name' returned for
%  each element in the object array.
            nout = min(numel(obj),narg);
            out = cell(nargout, nout);
            for co = 1:numel(obj)
                [out{:,co}] = feval(method_name, obj(co), varargin{:});
            end
        end

        function varargout = run_method(obj, method_name, varargin)
%   Run a method for all objects in object array obj
%
%   [...] = obj.run_method(method_name, ...) runs the method
%   specified by method_name for each object in the object array obj. The
%   input arguments should have the same size as the object array. Each
%   input element is passed to the function and the outputs returned are
%   cells holding the output of the function. The cells have the same size
%   as the object array. If all the outputs of the function are scalars, an
%   array is produced instead of a cell
%
%   [...] = obj.run_method(method_name, ..., 'NoExpand', ...) Anything
%   passed after the 'NoExpand' word is not expaned but passed as is to
%   each call of the function.
%
%   See also:
%   ArraySupport, assign_paramater
                siz_out = size(obj);
                argout = cell(numel(obj),nargout);
                if ~isempty(varargin)
                    % handle inputs that should not be extended
                    fnoexp = find(cellfun(@(x) ischar(x) &&...
                        strcmp(x, 'NoExpand'), varargin),...
                        1, 'first');
                    no_exp_in = varargin(fnoexp+1:end);
                    varargin(fnoexp:end) = [];

                    % handle scalar inputs
                    fscalar = cellfun(@(x) isequal(obj.get_size(x),[1,1]),varargin);
                    scal_in =  varargin(fscalar);
                    varargin(fscalar)=[];

                    % handle other input
                    if ~isempty(varargin)
                        siz_in = cellfun(@obj.get_size, varargin, ...
                            "UniformOutput",false);
                        assert(isempty(siz_in) || isequal(size(obj),...
                            siz_in{:}),...
                            'helper.ArraySupport:NonMatchingInputSize',...
                            ['Size of non-scalar inputs must match', ...
                            ' size of object array'])
                        varargin = cellfun(@(x) reshape(x,[],1),...
                            varargin, ...
                            'UniformOutput',false);
                        not_cell = ~cellfun(@iscell,varargin);
                        varargin(not_cell) = cellfun(@num2cell, ...
                            varargin(not_cell),'UniformOutput',false);
                        varargin = [varargin{:}];
                    else
                        varargin = cell.empty(numel(obj),0);
                    end
                else 
                    varargin = cell.empty(numel(obj),0);
                    scal_in = {};
                    no_exp_in = {};
                end
                for co = 1:numel(obj)
                    [argout{co,:}] = feval(method_name, obj(co),...
                        varargin{co,:}, scal_in{:}, no_exp_in{:});
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
        function val = get_size(x)
            if ischar(x)
                val = [1, 1];
            else
                val = size(x);
            end
        end
    end
end
