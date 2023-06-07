classdef EnsembleFilter < Filter
% Class to filter out ensembles
%   
%   obj = EnsembleFilter() constructs a default EnsembleFilter object
%
%   obj = EnsembleFilter(ADCP) constructs an EnsembleFilter object with
%   bad_ensembles set to false for all ensembles in ADCP
%
%   obj = EnsembleFilter(is_bad) constucts and EnsembleFilter object with
%   bad_ensembles set to the values given in is_bad. is_bad is a logical
%   row vector that indicates the ensembles that should be filtered out.
%
%   EnsembleFilter properties:
%   bad_ensembles - marks the bad ensembles
%
%   see also: Filter, cross_section_selector
    properties
        bad_ensembles (1,:) logical = logical.empty(1,0)
    end
    methods
        function obj = EnsembleFilter(varargin)
            obj = obj@Filter(varargin{:});
            varargin = obj(1).unprocessed_construction_inputs;
            for count_arg = 1 : numel(varargin)
                current_arg=varargin{count_arg};
                if isa(current_arg,'ADCP')
                    if ~isscalar(obj)
                        for co = 1:numel(obj)
                            obj(co).bad_ensembles=false(1,current_arg(co).nensembles);
                        end
                    else
                        obj.bad_ensembles=false(1,sum([current_arg.nensembles]));
                    end
                elseif islogical(current_arg)
                    obj.bad_ensembles=current_arg;
                end
            end
            [obj.description]=deal('Ensemble filter');
        end
    end
    methods(Access=protected)
        function bad=bad_int(obj,adcp)
            if isempty(obj.bad_ensembles)
                obj.bad_ensembles=false([1,sum([adcp.nensembles])]);
            end
           assert(numel(obj.bad_ensembles)==sum([adcp.nensembles]), 'The number of elements in the bad_ensembles property should match the number of ensembles in the ADCP object') 
           bad=repmat(obj.bad_ensembles,[max([adcp.ncells]),1,max([adcp.nbeams])]);
        end
    end
end