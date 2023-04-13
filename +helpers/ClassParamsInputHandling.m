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
            props_class = cellfun(@(x) class(obj(1).(x)),propnames,...
                'UniformOutput',false);
            [props_class,idx] = unique(props_class);
            propnames = propnames(idx);

            %%% find parameter inputs and assign them

            % find strign or character inputs that match assignable
            % property names
            fpars = find(cellfun(@(x) (isstring(x) || ischar(x)) &&...
                any(strcmp(x,propnames)),varargin));

            % check paramter names are not consecutive
            fpars(diff(fpars)<2) = [];

            % check parameter name is not last input
            fpars(fpars == numel(varargin))=[];

            fpars_success_pars = true(size(fpars));
            for ia = 1:numel(fpars)
                try
                    obj.assign_property(varargin{fpars(ia)},varargin{fpars(ia)+1});
                catch % if assignment is unsuccesfull, remove it from
                    % successfull assignments
                    fpars_success_pars(ia)=false;
                end
            end
            fpars_success_pars = fpars(fpars_success_pars);
            fpars_success_pars = [fpars_success_pars fpars_success_pars+1];

            %%% find class based inputs and assign them
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
                    obj.assign_property(propnames{fmatch},varargin{f_unproc(cp)})
                else
                    fpars_success_class(cp) = false;
                end
            end
            fpars_success_class = f_unproc(fpars_success_class);

            %%% get index of unhandled inputs
            unhandled = setdiff(1:numel(varargin),...
                [fpars_success_pars fpars_success_class]);

            %%% Generate warning or error if needed
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
            if nargout > 0
                out = unhandled;
            end
        end
    end
end
