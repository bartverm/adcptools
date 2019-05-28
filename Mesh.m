classdef Mesh < XSection
    properties
       cells(:,1) MeshCell = MeshCell.empty(0,1);
    end
    properties(Dependent, SetAccess=protected, GetAccess=public)
       n_cells;
    end
    properties(Dependent)
       coordinates(2,:,:) double {mustBeFinite} 
    end
    methods
        function val=get.n_cells(obj)
            val=size(obj.cells,1);
        end
        function val=get.coordinates(obj)
            val=cat(3,obj.cells.coordinates);
        end
        function set.coordinates(obj,val)
            assert(size(val,3)==obj.n_cells,'Size of third dimension must match number of cells in mesh')
            cors=squeeze(mat2cell(val,size(val,1),size(val,2),ones(size(val,3),1)));
            [obj.cells.coordinates]=deal(cors{:});
        end
        function varargout=plot(obj,varargin)
            hbed=plot@XSection(obj,'k','Linewidth',2);
            hcells=plot(obj.cells);
            if nargout > 0
                varargout{1}=hbed;
                varargout{2}=hcells;
            end
        end
    end
end