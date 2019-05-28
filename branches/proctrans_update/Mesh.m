classdef Mesh < XSection
    properties
       cells(:,1) MeshCell = MeshCell.empty(0,1);
    end
    properties(Dependent, SetAccess=protected, GetAccess=public)
       n_cells=size(cells,1);
    end
end