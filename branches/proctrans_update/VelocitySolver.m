classdef VelocitySolver < handle
    properties
        adcp (1,1) VMADCP
        mesh (1,1) Mesh = SigmaZetaMesh;
        bathy (1,1) Bathymetry = BathymetryScatteredPoints;
        xs (1,1) XSection
        filter (:,1) Filter
    end
    methods
        function obj=VelocitySolver(varargin)
            has_vmadcp=false;
            has_mesh=false;
            has_bathy=false;
            has_xs=false;
            for cnt_arg=1:nargin
                cur_arg=varargin{cnt_arg};
                if isa(cur_arg,'VMADCP')
                    has_vmadcp=true;
                    obj.adcp=cur_arg;
                elseif isa(cur_arg, 'Mesh')
                    has_mesh=true;
                    obj.mesh=cur_arg;
                elseif isa(cur_arg, 'Bathymetry')
                    has_bathy=true;
                    obj.bathy=cur_arg;
                elseif isa(cur_arg,'Filter')
                    obj.filter=cur_arg;
                elseif isa(cur_arg,'XSection')
                    has_xs=true;
                    obj.xs=cur_arg;
                end
            end
            if ~has_mesh
                error('You must provide a Mesh object upon construction of a VelocitySolver object')
            end
            if ~has_vmadcp
                error('You must provide a VMADCP object upon construction of a VelocitySolver object')
            end
            if ~has_bathy
                obj.bathy=BathymetryScatteredPoints(obj.adcp);
            end
            if ~has_xs
                obj.xs=XSection(obj.adcp);
            end
        end
    end
    methods (Abstract)
        vel=get_velocity(obj)
    end
end