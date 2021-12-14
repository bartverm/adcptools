classdef SigmaZetaMeshGenerator < handle
% Base class to produce SigmaZetaMeshes
%
%   Subclasses should implement the get_mesh method.
%
%   SigmaZetaMeshGenerator methods:
%   get_mesh - returns the generated mesh
%
%   see also: Mesh, SigmaZetaMeshFromVMADCP
    methods(Abstract)
        mesh=get_mesh(obj)
    end
end