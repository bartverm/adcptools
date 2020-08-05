classdef SigmaZetaMeshGenerator < handle
    methods(Abstract)
        mesh=get_mesh(obj)
    end
end