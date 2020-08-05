classdef EnsembleFilter < Filter
    % Class to filter out ensembles
    %
    properties
        bad_ensembles (1,:) logical
    end
    methods
        function obj=EnsembleFilter(varargin)
            for count_arg=1:nargin
                current_arg=varargin{count_arg};
                if isa(current_arg,'ADCP')
                    obj.bad_ensembles=false(1,current_arg.nensembles);
                elseif islogical(current_arg)
                    obj.bad_ensembles=current_arg;
                end
            end
            obj.description='Ensemble filter';
        end
    end
    methods(Access=protected)
        function bad=bad_int(obj,adcp)
           assert(numel(obj.bad_ensembles)==adcp.nensembles, 'The number of elements in the bad_ensembles property should match the number of ensembles in the ADCP object') 
           bad=repmat(obj.bad_ensembles,[max(adcp.ncells),1,max(adcp.nbeams)]);
        end
    end
end