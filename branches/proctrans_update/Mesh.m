classdef Mesh < handle
% Base class for Meshes for ADCP data processing
%
%   Subclasses need to implement the get_ncells, index, and plot methods
%
%   Mesh properties (read only):
%   ncells - number of cells
%
%   Mesh methods:
%   index - returns mesh cell indices of cells given points fall into
%   plot - plot the mesh optionally coloring with a given variable
%
%   see also: SigmaZetaMesh, VelocitySolver
    properties (Dependent)
        ncells
    end
    methods
        function val=get.ncells(obj)
            val=obj.get_ncells();
        end
    end
    methods (Abstract)
        index(obj,n, sigma)
        plot(obj,var)
    end
    methods(Access=protected, Abstract)
        get_ncells(obj)
    end
end