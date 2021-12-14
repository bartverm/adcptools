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
    methods
        function             obj=LocationBasedVelocitySolver(varargin)
            obj=obj@VelocitySolver(varargin{:});
            for cnt_arg=1:nargin
                cur_arg=varargin{cnt_arg};
                if isa(cur_arg, 'VMADCP') || isa(cur_arg, 'Mesh') || isa(cur_arg,'Bathymetry') || isa(cur_arg, 'Filter') || isa(cur_arg, 'XSection')
                    continue
                else
                    warning(['Unhandled input of type: ', class(cur_arg), ' on construction of LocationBasedVelocitySolver object'])
                end
            end
        end
    end
    methods(Access=protected)
        function [vpos, vdat, xform] = get_solver_input(obj)
            % Get velocity position and compute sigma coordinates
            vpos=obj.adcp.depth_cell_position; % velocity positions
            
            % get velocity data
            vdat = obj.adcp.water_velocity(CoordinateSystem.Beam); % get velocity data

            % get transformation matrix
            xform = obj.adcp.xform(CoordinateSystem.Beam,CoordinateSystem.Earth); % get Earth to Beam transformation matrix
            xform(:,:,:,4)=[]; % remove Error velocity to beam transformation
        end

    end

end