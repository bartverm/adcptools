classdef LocationBasedVelocitySolver < VelocitySolver
    % Implements a time based velocity solver
    %
    %   This class solves the water velocity on the given mesh, by 
    %   combining beam velocities based on the location they were measured.
    %   Beam velocities measure in the same cell are used to fit a least
    %   squares velocity vector, taking into account the tilting of the
    %   ship
    %
    %   obj=LocationBasedVelocitySolver(...) allows to set properties upon 
    %   construction.
    %   Depending on the class objects are assigned to a specific property:
    %   VMADCP - assigned to the adcp property
    %   Mesh - assigned to the mesh property
    %   Bathymetry - assigned to the bathymetry property
    %   XSection - assigned to the xs property
    %   Filter - assigned to the filter property
    %
    %  LocationBasedVelocitySolver properties:
    %   adcp - VMADCP object with the adcp data
    %   mesh - Mesh object defining a mesh
    %   bathy - Defines the bathymetry
    %   xs - Defines the cross-section
    %   filter - Selects ADCP data to be excluded from the computations
    %
    %   LocationBasedVelocitySolver methods:
    %   get_parameters - get velocity model parameters on mesh
    %   get_velocity - get velocity on mesh
    %   rotate_to_xs - rotates velocity to cross-section direction
    %
    %   see also: VMADCP, Mesh, Bathymetry, XSection, Filter,
    %   VelocitySolver
    methods(Access=protected)
        function [vpos, vdat, xform, time, wl] = get_solver_input(obj)
            [vpos, ~, ~, time, wl] = get_solver_input@ADCPDataSolver(obj);

            % get velocity data
            vdat = obj.adcp.cat_property('water_velocity',...
                CoordinateSystem.Beam); % get velocity data in Beam coords

            % get transformation matrix
            xform = obj.adcp.cat_property('xform',...
                CoordinateSystem.Beam,... destination: Beam
                CoordinateSystem.Earth,...% source: Earth
                'Geometry', true); % Force use of tilts
            xform(:,:,:,4)=[]; % remove Error velocity to beam transformation
            
            % filter and vectorize
            [vdat, xform] = obj.filter_and_vectorize(vdat, xform);
        end
    end

end
