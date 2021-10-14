classdef VelocitySolver < handle
    % Base class to solve ADCP repeat transect velocity data
    %
    %   Subclasses should implement the get_velocity method. This method
    %   produces an array of velocities in each of the mesh cells.
    %
    %   obj=VelocitySolver(...) allows to set properties upon construction.
    %   Depending on the class objects are assigned to a specific property:
    %   VMADCP - assigned to the adcp property
    %   Mesh - assigned to the mesh property
    %   Bathymetry - assigned to the bathymetry property
    %   XSection - assigned to the xs property
    %   Filter - assigned to the filter property
    %
    %  VelocitySolver properties:
    %   adcp - VMADCP object with the adcp data
    %   mesh - Mesh object defining a mesh
    %   bathy - Defines the bathymetry
    %   xs - Defines the cross-section
    %   filter - Selects ADCP data to be excluded from the computations
    %
    %   VelocitySolver methods:
    %   get_velocity - returns velocities for each mesh cell (Abstract)
    %   get_velocity_xs - rotates velocity to cross-section direction
    %
    %   see also: VMADCP, Mesh, Bathymetry, XSection, Filter,
    %   TimeBasedVelocitySolver
    properties
        % VelocitySolver/adcp
        %
        %   Scalar VMADCP object holding the adcp data to compute the velocity
        %
        %   see also: VelocitySolver, VMADCP
        adcp (1,1) VMADCP
        
        % VelocitySolver/mesh
        %
        %   Scalar Mesh object defining the mesh on which the data is to be solved.
        %   Default value is SigmaZetaMesh
        %
        %   see also: VelocitySolver, Mesh
        mesh (1,1) Mesh = SigmaZetaMesh;
        
        % VelocitySolver/bathy
        %
        %   Scalar bathymetry object that defines the location of the bed. Default
        %   is BathymetryFromScatteredPoints(adcp), which constructs a bathymetry
        %   based on the VMADCP data in adcp.
        %
        %   see also: VelocitySolver, BathymetryScatteredPoints
        bathy (1,1) Bathymetry = BathymetryScatteredPoints;
        
        % VelocitySolver/xs
        %
        %   Scalar XSection object defining the origin and direction of the
        %   cross-section. Default value is XSection(adcp) which construct a
        %   cross-section based on the track of the VMADCP data in adcp
        %
        %   see also: VelocitySolver, XSection
        xs (1,1) XSection
        
        % VelocitySolver/filter
        %
        %   Vector with filters to exclude data from the velocity processing.
        %
        %   see also: VelocitySolver, Filter
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
        [vel, velstd]=get_velocity(obj)
    end
    methods
        function vel=get_velocity_xs(obj)
        % Rotates velocity to direction of cross-section
        %
        %   vel=get_velocity_xs(obj) Rotates the velocity obtained
        %   with get_velocity to the direction of the cross-section.
        %
        % See also: VelocitySolver, get_velocity
            orig_vel=obj.get_velocity();
            [~,~,vels,veln]=obj.xs.xy2sn(obj.mesh.x_middle(obj.mesh.col_to_cell),obj.mesh.y_middle(obj.mesh.col_to_cell),orig_vel(:,1),orig_vel(:,2));
            vel=[vels,veln,orig_vel(:,3)];
        end
    end
end