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
    %   Filter - assigned to the ensemble_filter property
    %
    %  VelocitySolver properties:
    %   adcp - VMADCP object with the adcp data
    %   mesh - Mesh object defining a mesh
    %   bathy - Defines the bathymetry
    %   xs - Defines the cross-section
    %   ensemble_filter - defines repeat transects to be processed
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
        
        % VelocitySolver/mesh mesh on which velocity is solved
        %
        %   Mesh object defining the mesh on which the data is to be 
        %   solved. If an array is passed, each mesh is used for a
        %   different repeat transect, thus the number of elements must
        %   match the number of elements in the ensemble_filter property.
        %   
        %   Default value is SigmaZetaMesh
        %
        %   see also: VelocitySolver, Mesh
        mesh (1,:) Mesh = SigmaZetaMesh;
        
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
        
        % VelocitySolver/ensemble_filter
        %
        %   Ensemble filters defining repeat crossings to be processed
        %
        %   see also: VelocitySolver, EnsembleFilter
        ensemble_filter (1,:) EnsembleFilter
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
                    obj.ensemble_filter=cur_arg;
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
% returns velocity solved on mesh cells
%
%   vel=obj.get_velocity() returns the velocity as an Nx3xM matrix with N 
%       being the number of cells in the mesh, 3 is the number of
%       components solved, and M is the number of repeat crossings, which
%       corresponds to the number of elements in EnsembleFilter
%   [vel, velstd]=obj.get_velocity() also returns the standard deviation in
%       the solved velocity data. Dimensions are the same as the velocity
%       output
%
%   see also: VelocitySolver, get_velocity_xs
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
            idx_mesh=obj.make_indices();
            orig_vel=obj.get_velocity();
            vel=cell(size(orig_vel));
            for crp=1:numel(orig_vel)
                cur_mesh=obj.mesh(idx_mesh(crp));
                [~,~,vels,veln]=obj.xs.xy2sn(cur_mesh.x_middle(cur_mesh.col_to_cell),cur_mesh.y_middle(cur_mesh.col_to_cell),orig_vel{crp}(:,1),orig_vel{crp}(:,2));
                vel{crp}=[vels,veln,orig_vel{crp}(:,3)];
            end
        end
    end
    methods(Access=protected)
        function [idx_mesh,idx_ef]=make_indices(obj)
            nrp=max(numel(obj.ensemble_filter),numel(obj.mesh));
            idx_mesh=ones(nrp,1);
            idx_ef=ones(nrp,1);
            if ~isscalar(obj.mesh)
                idx_mesh=cumsum(idx_mesh);
                assert(numel(obj.mesh)==nrp,'number of elements of mesh and ensemble_filter properties should match')
            end
            if ~isscalar(obj.ensemble_filter)
                idx_ef=cumsum(idx_ef);
                assert(numel(obj.ensemble_filter)==nrp,'number of elements of mesh and ensemble_filter properties should match')
            end
        end
    end
end