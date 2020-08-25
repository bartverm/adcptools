classdef TimeBasedVelocitySolver < VelocitySolver
    % Implements a time based velocity solver
    %
    %   This class solves the water velocity on the given mesh, by first
    %   computing the velocity from the ADCP data in a conventional way, i.e.
    %   by combining beam velocities based on the time they were measured.
    %   Subsequently, velocities measured in the same cell are averaged to
    %   produce the output velocity. The positions of the velocities are
    %   matched in sigma coordinates.
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
    %   get_velocity - returns velocities for each mesh cell
    %
    %   see also: VMADCP, Mesh, Bathymetry, XSection, Filter,
    %   VelocitySolver
    methods
        function obj=TimeBasedVelocitySolver(varargin)
            obj=obj@VelocitySolver(varargin{:});
            for cnt_arg=1:nargin
                cur_arg=varargin{cnt_arg};
                if isa(cur_arg, 'VMADCP') || isa(cur_arg, 'Mesh') || isa(cur_arg,'Bathymetry') || isa(cur_arg, 'Filter') || isa(cur_arg, 'XSection')
                    continue
                else
                    warning(['Unhandled input of type: ', class(cur_arg), ' on construction of TimeBasedVelocitySolver object'])
                end
            end
        end
        function vel=get_velocity(obj)
            % Solves velocity combining beam velocities based on time
            %
            %   vel=get_velocity(obj) computes the velocity in the mesh cells by using
            %   the velocity computed by combining the beam velocities measured at the
            %   same time (this is done by the ADCP class). Velocities that where
            %   measured within a mesh cell (this is determined by the Mesh class) are
            %   averaged.
            %
            % see also: TimeBasedVelocitySolver, Mesh, VMADCP
            vpos=obj.adcp.depth_cell_position;
            filters=[obj.filter; obj.adcp.filters];
            vpos(filters.bad(obj.adcp))=nan;
            vpos=nanmean(vpos,3);
            [~, n_pos]=obj.xs.xy2sn(vpos(:,:,:,1),vpos(:,:,:,2));
            zb_pos=obj.bathy.get_bed_elev(vpos(:,:,:,1), vpos(:,:,:,2));
            sig_pos=1-vpos(:,:,:,3)./zb_pos;
            idx=repmat(obj.mesh.index(n_pos, sig_pos), 1, 1, 3);
            fcomp=cumsum(ones(size(idx)),3);
            vel=nan(obj.mesh.ncells, 3);
            vel_data = obj.adcp.water_velocity(CoordinateSystem.Earth);
            for cc = 1:obj.mesh.ncells
                for comp=1:3
                    vel(cc,comp) = nanmean(vel_data(idx==cc & fcomp==comp));
                end
            end
        end
    end
end