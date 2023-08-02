classdef ArraySupport < handle
% Class to help onstruct object arrays based on array input to constructor
%     
%   obj = ArraySupport(...) constructs the object array obj. The size of 
%   the object array is determined based on the inputs, and the inputs are
%   assigned to the objects in the array. Which inputs are considered to
%   determine the size and how they are assigned to the array depends on
%   the type of inputs:
%   name = value inputs are not included in determining the size of the
%       output array. The given value is assigned as is to each of the
%       objects in the array.
%   'ParamName', paramValue input pairs are considered in determining the
%       size only if 'ParamName' is an assignable property of the class and
%       the pair is not given after the 'NoExpand' input argument (see
%       below). The paramValue is assigned to object property by assigning
%       each element of paramValue to the property with name 'ParamName' of
%       the element in the object array.
%   handleClasses are used to determine the output array size only if they
%       are assignable based on their class, i.e. there is only one
%       property with that particular handle class. The elements in the
%       input array are assigned to the property with corresponding class
%       of the corresponding element in the output array
%   cell inputs are used to determine the size of the output array only if
%       all the values held by the cell are of the same assignable handle 
%       class.
%   All inputs that determine the size of the output array must all be
%   scalar or have matching size.
%       
%   obj = ArraySupport('NoExpand',...) any input given after this argument
%   is not used to determine size of constructed array. All name=value
%   pairs cannot be given before the 'NoExpand' argument and therefor will
%   never be expanded.
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
    properties(Access = protected, Transient)
        unprocessed_construction_inputs (1,:) cell = {};
        unassigned_properties (1,:) cell = {};
    end
    methods
        function obj = ArraySupport(varargin)
            if isempty(varargin)
                return
            end

            % inputs to use to determine size
            expand_inputs = false(size(varargin)); % variables to expand
            processed_inputs = false(size(varargin)); % used inputs
            propname = cell(size(varargin)); % property name to assign 
                       
            % keep value of all assignable 'ParamName', paramValue pairs
            f_unprocessed_inputs = find(~processed_inputs);
            f_pnampval = obj.find_paramname_paramvalue(...
                varargin{~ processed_inputs});
            f_pnampval = f_unprocessed_inputs(f_pnampval);
            propname(f_pnampval+1) = varargin(f_pnampval);
            processed_inputs(f_pnampval) = true;
            processed_inputs(f_pnampval + 1) = true;
            expand_inputs(f_pnampval + 1) = true;

            % keep value of all assignable name=value pairs
            f_unprocessed_inputs = find(~processed_inputs);
            f_namval = obj.find_name_value(varargin{~processed_inputs});
            f_namval = f_unprocessed_inputs(f_namval);
            processed_inputs(f_namval) = true;
            processed_inputs(f_namval + 1) = true;
            propname(f_namval + 1) = varargin(f_namval);
            propname(f_namval + 1) = cellfun(@char,...
                propname(f_namval + 1), UniformOutput=false);

            % cell inputs
            f_unprocessed_inputs = find(~processed_inputs);
            [f_cell_inputs, assignable, cell_names] = ...
                obj.find_cell(varargin{~processed_inputs});
            f_cell_inputs = f_unprocessed_inputs(f_cell_inputs);
            expand_inputs(f_cell_inputs) = true;
            processed_inputs(f_cell_inputs(assignable)) = true;
            propname(f_cell_inputs(assignable)) = cell_names;

            % handle class inputs
            f_unprocessed_inputs = find(~processed_inputs);
            [f_handle_inputs, assignable, handle_names] =...
                obj.find_handle(varargin{~processed_inputs});
            f_handle_inputs = f_unprocessed_inputs(f_handle_inputs);
            expand_inputs(f_handle_inputs) = true;
            processed_inputs(f_handle_inputs(assignable)) = true;
            propname(f_handle_inputs(assignable)) = handle_names(assignable);

            % remove all inputs given after 'NoExpand'
            f_expand = find(strcmp('NoExpand',varargin),1,'first');
            if ~isempty(f_expand) 
                expand_inputs(f_expand:end) = false;
            end

            % don't expand scalar inputs
            f_scalar = cellfun(@isscalar, varargin);
            expand_inputs(f_scalar) = false;
    
            %%% create object array
            if any(expand_inputs)
                % get size of expanded variables
                siz_nscal = cellfun(@size,varargin(expand_inputs), ...
                    'UniformOutput', false);
                % make sure their sizes match
                assert(isscalar(siz_nscal) || isequal(siz_nscal{:}),...
                    'Solver:NonMatchingInputSize',...
                    'Size of all non-scalar input should match')
                siz_obj = num2cell(siz_nscal{1});
                obj(siz_obj{:}) = feval(class(obj(1)));
                for co = 2:numel(obj)-1 % needed for heterogeneous classes
                    obj(co) = feval(class(obj(1)));
                end
            end

            %%% assign inputs
            f_remove = cellfun(@isempty, propname);
            propname(f_remove) = [];
            assign_vars = varargin;
            assign_vars(f_remove)=[];
            expand_inputs(f_remove) = [];
            exp_states = {'NoExpand',[]};
            expand_inputs = exp_states(expand_inputs + 1);
            for ci = 1:numel(propname)
                obj.assign_property(propname{ci}, assign_vars{ci},...
                    expand_inputs{ci})
            end
            allprops = obj.find_assignable_properties();
            obj(1).unassigned_properties = setdiff(allprops, propname);
            obj(1).unprocessed_construction_inputs = ...
                varargin(~ processed_inputs);

            obj.init_handle_properties();
        end
    end
    methods(Access = public, Sealed)
        function init_handle_properties(obj)
            if isscalar(obj)
                return
            end
            propnames = obj(1).unassigned_properties;
            f_handle = cellfun(@(x) isa(obj(1).(x), 'handle'), propnames);
            propnames = propnames(f_handle);
            for co=2:numel(obj)
                for cp=1:numel(propnames)
                    obj(co).(propnames{cp}) =...
                        feval(class(obj(1).(propnames{cp})));
                end
            end
        end
        function propnames = find_assignable_properties(obj, access_class)
            if nargin < 2 % if access object is given
                access_class = obj; % if no access class is given check for 
                    % itself
            end

            mc = metaclass(obj);
            propnames = {mc.PropertyList.Name};
            f_props =...
                ... access is one stated in access_list
                (cellfun(@(x) ischar(x) &&...
                ismember(x, {'public'}),...
                {mc.PropertyList.SetAccess}) | ... 
                ... obj is of class that has access
                cellfun(@(x) iscell(x) &&...
                any(cellfun(@(y) isa(access_class,y.Name), x)),...
                {mc.PropertyList.SetAccess})) & ...
                ... properties are not dependent
                ~[mc.PropertyList.Dependent] & ...
                ... properties are not constant
                ~[mc.PropertyList.Constant] ;
            propnames = propnames(f_props);
        end
        function [propnames, props_class] =...
            find_assignable_handle_properties(obj)
            propnames = obj.find_assignable_properties();
            if isempty(propnames)
                props_class ={};
                return
            end
            mc = metaclass(obj);
            plist = mc.PropertyList;
            allprops = {plist.Name};
            [~, f_assignable] = intersect(allprops,propnames,'stable');
            plist = plist(sort(f_assignable));
            val = [plist.Validation];
            val_class = [val.Class];
            props_class = {val_class.Name};
            [props_class,idx] = unique(props_class,'stable');
            propnames = propnames(idx);
        end
        function fpars = find_paramname_paramvalue(obj,varargin)
            propnames = obj.find_assignable_properties;

            % find string or character inputs that match assignable
            % property names
            fpars = find(cellfun(@(x) ischar(x) &&...
                any(strcmp(x,propnames)),varargin));

            % check paramter names are not consecutive
            fpars(diff(fpars)<2) = [];

            % check parameter name is not last input
            fpars(fpars == numel(varargin))=[];
        end

        function fpars = find_name_value(obj,varargin)
            propnames = obj.find_assignable_properties;

            % find string or character inputs that match assignable
            % property names
            fpars = find(cellfun(@(x) isstring(x) &&...
                any(strcmp(x,propnames)),varargin));

            % check paramter names are not consecutive
            fpars(diff(fpars)<2) = [];

            % check parameter name is not last input
            fpars(fpars == numel(varargin))=[];
        end
        
        function [f_handle, assignable, prop_names] = find_handle(obj,varargin)
            [hpropnames, hpropclass] = obj.find_assignable_handle_properties();
            f_handle = find(cellfun(@(x) isa(x,'handle'), varargin));
            handle_inputs = varargin(f_handle);
            assignable = false(size(handle_inputs));
            prop_names = cell(size(handle_inputs));
            for ci = 1:numel(handle_inputs)
                fprop=cellfun(@(x) isa(handle_inputs{ci},x), hpropclass);
                if any(fprop)
                    assignable(ci) = true;
                    fprop = find(fprop,1,'first');
                    prop_names(ci) = hpropnames(fprop);
                end
            end
        end
        function [f_cell, assignable, prop_names] = find_cell(obj, varargin)
            [hpropnames, hpropclass] = obj.find_assignable_handle_properties();
            f_cell = find(cellfun(@iscell, varargin));
            cell_inputs = varargin(f_cell);
            prop_names = cell(size(cell_inputs));
            cell_is_assignable = false(size(cell_inputs));
            for cc = 1:numel(f_cell)
                cur_cell = cell_inputs{cc};
                if ~all(cellfun(@(x) isa(x,'handle'), cur_cell))
                    continue
                end
                cell_classes = cellfun(@class, cur_cell,...
                    'UniformOutput',false);
                if ~isequal(cell_classes{:})
                    continue
                end
                % cur_class = cell_classes{1};
                fmatch = cellfun(@(x) isa(cur_cell{1}, x), hpropclass);
                if any(fmatch)
                    prop_names(cc) = hpropnames(fmatch);
                    cell_is_assignable(cc) = true;
                end
            end
            assignable = cell_is_assignable;
            prop_names = prop_names(cell_is_assignable);
        end
        function warn_unprocessed(obj)
            unhandled = obj(1).unprocessed_construction_inputs;
            if ~isempty(unhandled)
                if numel(unhandled) == 1
                    inp_num_str = sprintf('%d',unhandled);
                else
                    inp_num_str = sprintf('%d, ',unhandled(1:end-1));
                    inp_num_str(end-1:end)=[];
                    inp_num_str = [inp_num_str ' and ' sprintf('%d',unhandled(end))];
                end
                msg_id = 'ClassParamsInputHandling:unhadled_inputs';
                msg = ['Could not handle inputs number ', inp_num_str];
                warning(msg_id, msg)
            end
        end

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
