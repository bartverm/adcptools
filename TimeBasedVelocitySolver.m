classdef TimeBasedVelocitySolver < VelocitySolver
    % Implements a time based velocity solver
    %
    %   This class solves the water velocity on the given mesh, by first
    %   computing the velocity from the ADCP data in a conventional way, i.e.
    %   by combining beam velocities based on the time they were measured.
    %   Subsequently, velocities measured in the same cell are averaged to
    %   produce the output velocity. The positions of the velocities are
    %   matched in sigma coordinates. If a velocity model is given, the
    %   model parameters will be fit based on the cartesian component of
    %   velocity
    %
    %   obj=TimeBasedVelocitySolver(...) allows to set properties upon construction.
    %   Depending on the class objects are assigned to a specific property:
    %   VMADCP - assigned to the adcp property
    %   Mesh - assigned to the mesh property
    %   Bathymetry - assigned to the bathymetry property
    %   XSection - assigned to the xs property
    %   Filter - assigned to the filter property
    %
    %  TimeBasedVelocitySolver properties:
    %   adcp - VMADCP object with the adcp data
    %   mesh - Mesh object defining a mesh
    %   bathy - Defines the bathymetry
    %   xs - Defines the cross-section
    %   filter - Selects ADCP data to be excluded from the computations
    %
    %   TimeBasedVelocitySolver methods:
    %   get_parameters - get velocity model parameters on mesh
    %   get_velocity - get velocity on mesh
    %   rotate_to_xs - rotates velocity to cross-section direction
    %
    %   see also: VMADCP, Mesh, Bathymetry, XSection, Filter,
    %   VelocitySolver
    methods(Access=protected)
        function [vpos, vdat, xform, time, wl] = get_solver_input(obj)
            [vpos, ~, ~, time, wl] = get_solver_input@ADCPDataSolver(obj);

            % Get velocity position and compute sigma coordinates
            vpos = mean(vpos, 3, 'omitnan'); % average position of four beams
            vpos = repmat(vpos, [1, 1, 4, 1]); % replicate average position

            % get velocity data
            vdat = obj.adcp.cat_property('water_velocity',...
                CoordinateSystem.Earth); % get velocity data

            % get transformation matrix
            xform = obj.adcp.cat_property('xform',...
                CoordinateSystem.Earth,... destination: Earth
                CoordinateSystem.Earth); % source: Earth
            xform(:,:,:,4)=[]; % remove Error velocity

            % filter and vectorize
            [vdat, xform] = obj.filter_and_vectorize(vdat, xform);
        end
    end
end