classdef Filter <  handle & ...
        matlab.mixin.Heterogeneous & ...
        helpers.ArraySupport
% Generic class to implement filters for ADCP objects
%
%   This class implements a dummy filter, i.e. it does not filter anything.
%   To implement a usefull filter you can subclass this class and implement
%   the bad_int function. 
%   This class inherits from mixin.Heterogenous, which allows to create
%   arrays with subclasses of Filter. This allows to combine several
%   different filters when processing ADCP data. Calling the 'bad' method
%   on the object array, will return the cells 
%
%   Filter methods (Sealed):
%   bad - returns logical array marking bad cells
%   all_cells_bad - returns logical array marking ens with all bad cells
%   any_cells_bad - returns logical array marking ens with any bad cells

    properties(SetAccess=protected)
        % A brief description of the filter
        description(1,:) char='Dummy filter';
    end
    methods
        function obj = Filter(varargin)
            obj = obj@helpers.ArraySupport(varargin{:});
        end
    end
    methods(Sealed)
        function bad=all_cells_bad(obj,adcp)
        % bad_ensemble(obj,ADCP) get ensembles in which all cells are bad.
        % For arrays of filters it will return ensembles that have all
        % cells marked in at least one of the filters
            bad=obj.expand_bad(@(obj,adcp) all(all(obj.bad(adcp),1),3), adcp);
        end

        function bad=any_cells_bad(obj,adcp)
        % bad_ensemble(obj,ADCP) get ensembles in which at least one bad
        % cell. For arrays of filters it will return ensembles that have 
        % at least one bad cell in at least one of the filters
            bad=obj.expand_bad(@(obj,adcp) any(any(obj.bad(adcp),1),3), adcp);
        end
        function bad=bad(obj,adcp)
        % bad(obj,ADCP) get bad cells for ADCP object. If obj is an
        % array of filters it will return the combined filter
            bad=obj.expand_bad(@(obj,adcp) obj.bad_int(adcp), adcp);
        end
        function plot(obj,adcp)
            figure
            nb=max(adcp.nbeams);
            no=numel(obj);
            axh=nan(nb*(no+1));
            ca=0;
            for co = 1:no
                filt=obj(co).bad(adcp);
                for cb=1:nb
                    ca=ca+1;
                    axh(ca)=subplot(no+1,nb,(co-1)*nb+cb);
                    imagesc(filt(:,:,cb))
                    title(obj(co).description)
                end
            end
            filt=obj.bad(adcp);
            for cb=1:nb
                ca=ca+1;
                axh(ca)=subplot(no+1,nb,no*nb+cb);
                imagesc(filt(:,:,cb))
                title('All filters')
            end
            linkaxes(axh,'xy')
        end
    end
    
    methods (Access=protected)
        function bad=bad_int(obj,adcp)
            if isscalar(obj)
                bad=false(max(adcp.ncells),adcp.nensembles, max(adcp.nbeams));
            else
                bad=obj(1).bad();
            end
        end
    end
    methods(Sealed, Access=private)
        function bad=expand_bad(obj, func, adcp) 
            % combines filter functions on object arrays
            if isempty(obj)
                bad=false;
                return
            end
            if isscalar(obj)
                bad=func(obj,adcp);
            else
                bad=obj(1).expand_bad(func,adcp);
                for co=1:numel(obj)
                    bad = bad | obj(co).expand_bad(func, adcp);
                end
            end
        end
    end
end