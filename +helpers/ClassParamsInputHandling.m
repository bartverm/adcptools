classdef ClassParamsInputHandling < helpers.ArraySupport
% class adding support for class and parameter based inputs
%
%  ClassParamsInputHandling methods:
%   parse_class_params_inputs - parses inputs
%
%   see also: helpers

    methods
        function out = parse_class_params_inputs(obj,varargin)
% Assign inputs to the class properties based on parameter name or class
%
%   obj.parse_class_params_inputs(varargin) assigns any 'parameterName', 
%   parameterValue pair to an objects' property whenever the 
%   'parameterName' matches an assignable object property. All inputs that 
%   do not fullfill above requirements and that are handle objects, are 
%   assigned to objects property when their class matches an assignable 
%   property's class. 
%
%   unhandled = obj.parse_class_params_inputs(varargin) returns the number
%   of the unhandles input parameters.
%
%   The method warns for unhandled inputs
%
%   see also: ClassParamsInputHandling

            %%% find assignable properties in class
            mc = metaclass(obj);
            propnames = {mc.PropertyList.Name};
            f_props = strcmp({mc.PropertyList.SetAccess},'public') & ...
                ~[mc.PropertyList.Dependent] & ...
                ~[mc.PropertyList.Constant] ;
            propnames = propnames(f_props);

            %%% handle NoExpand keyword
            f_no_expand = find(cellfun(@(x) isa(x,'char') &&...
                strcmp(x, 'NoExpand'), varargin), 1, 'first');
            no_expand = cell(numel(varargin),1);
            if ~isempty(f_no_expand)
                no_expand{f_no_expand:end} = 'NoExpand';
            end

            %%% find parameter inputs and assign them

            % find string or character inputs that match assignable
            % property names
            fpars = find(cellfun(@(x) (isstring(x) || ischar(x)) &&...
                any(strcmp(x,propnames)),varargin));

            % check paramter names are not consecutive
            fpars(diff(fpars)<2) = [];

            % check parameter name is not last input
            fpars(fpars == numel(varargin))=[];

            fpars_success_pars = true(size(fpars));
            prop_assigned=false(size(propnames));
            for ia = 1:numel(fpars)
                try
                    obj.assign_property(varargin{fpars(ia)},...
                        varargin{fpars(ia)+1},...
                        no_expand{fpars(ia)});
                catch err% if assignment is unsuccesfull, remove it from
                    % successfull assignments
                    fpars_success_pars(ia)=false;
                    warning(['Failed to assign property %s with error: ',...
                        err.message],varargin{fpars(ia)})
                end
                if fpars_success_pars(ia)
                    prop_assigned = true;
                end
            end
            fpars_success_pars = fpars(fpars_success_pars);
            fpars_success_pars = [fpars_success_pars fpars_success_pars+1];

            %%% find class based inputs and assign them
            propnames(prop_assigned) = []; % remove assigned properties from assignable ones

            % only assign properties based on class if there is only one
            % property with that particular class
            props_class = cellfun(@(x) class(obj(1).(x)),propnames,...
                'UniformOutput',false);
            [props_class,idx] = unique(props_class);
            propnames = propnames(idx);
             

            f_unproc = setdiff(1:numel(varargin), fpars_success_pars);
            fpars_success_class = true(size(f_unproc));
            % assign inputs based on class
            for cp = 1:numel(f_unproc)
                if ~isa(varargin{f_unproc(cp)},'handle')
                    fpars_success_class(cp) = false;
                    continue
                end
                fmatch = cellfun(@(x) isa(varargin{f_unproc(cp)},x),props_class);
                if any(fmatch)
                    obj.assign_property(propnames{fmatch},...
                        varargin{f_unproc(cp)},...
                        no_expand{f_unproc(cp)})
                else
                    fpars_success_class(cp) = false;
                end
            end
            fpars_success_class = f_unproc(fpars_success_class);

            %%% handle cell inputs if cell contains objects of class that 
            %%% is assignable
            f_unproc = setdiff(1:numel(varargin), ...
                [fpars_success_pars, fpars_success_class]);
            f_cell = f_unproc(cellfun(@iscell, varargin(f_unproc)));
            fpars_success_cell=true(size(f_unproc));
            for cc = 1:numel(f_cell)
                cur_cell = varargin{f_cell(cc)};
                if ~all(cellfun(@(x) isa(x,'handle'), cur_cell))
                    fpars_success_cell(cc) = false;
                    continue
                end
                cell_classes = cellfun(@class, cur_cell, 'UniformOutput',false);
                if ~isequal(cell_classes{:})
                    fpars_success_cell(cc) = false;
                    continue
                end
                cur_class = cell_classes{1};
                fmatch = strcmp(cur_class, props_class);
                if any(fmatch)
                    obj.assign_property(propnames{fmatch},...
                        cur_cell,...
                        no_expand{f_cell(cc)})
                else
                    fpars_success_cell(cc) = false;
                end
            end
            fpars_success_cell= f_unproc(fpars_success_cell);

            %%% get index of unhandled inputs
            unhandled = setdiff(1:numel(varargin),...
                [fpars_success_pars,...
                fpars_success_class,...
                fpars_success_cell]);

            %%% Generate warning or error if needed
            % if ~isempty(unhandled)
            %     if numel(unhandled) == 1
            %         inp_num_str = sprintf('%d',unhandled);
            %     else
            %         inp_num_str = sprintf('%d, ',unhandled(1:end-1));
            %         inp_num_str(end-1:end)=[];
            %         inp_num_str = [inp_num_str ' and ' sprintf('%d',unhandled(end))];
            %     end
            %     msg_id = 'ClassParamsInputHandling:unhadled_inputs';
            %     msg = ['Could not handle inputs number ', inp_num_str];
            %     warning(msg_id, msg)
            % end
            if nargout > 0
                out = unhandled;
            end
        end
    end
end
