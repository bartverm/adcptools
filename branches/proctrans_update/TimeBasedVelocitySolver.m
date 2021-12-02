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
        function [vel, velstd]=get_velocity(obj)
            % Solves velocity combining beam velocities based on time
            %
            %   vel=get_velocity(obj) computes the velocity in the mesh cells by using
            %   the velocity computed by combining the beam velocities measured at the
            %   same time (this is done by the ADCP class). Velocities that where
            %   measured within a mesh cell (this is determined by the Mesh class) are
            %   averaged.
            %
            %   [vel,std] = get_velocity(obj) also returns the standard
            %   deviation in the velocity
            %
            % see also: TimeBasedVelocitySolver, Mesh, VMADCP

            [idx_mesh,idx_ef]=obj.make_indices();

            % Sigma coordinates
            vpos=obj.adcp.depth_cell_position; % velocity positions
%             vpos(obj.adcp.filters.bad(obj.adcp))=nan; % filter out velocity (make this optional?)
            vpos=mean(vpos,3,'omitnan'); % compute average position of beams
            [~, n_pos]=obj.xs.xy2sn(vpos(:,:,:,1),vpos(:,:,:,2)); % get velocity position transverse to cross section
            zb_pos=obj.bathy.get_bed_elev(vpos(:,:,:,1), vpos(:,:,:,2)); % get bed elevation at velocity position. This could optionally be done with the bed detection at the beam.
            sig_pos=(vpos(:,:,:,3)-zb_pos)./(obj.adcp.water_level-zb_pos); % compute sigma coordinate of velocities

            % get velocity data
            vel_data = obj.adcp.water_velocity(CoordinateSystem.Earth); % get velocity data
            vel_data(:,:,4)=[]; % remove error velocity

            vel=cell(numel(idx_ef),1);
            velstd=cell(numel(idx_ef),1);
            for crp=1:numel(idx_ef)
                cur_n=n_pos(:,~obj.ensemble_filter(idx_ef(crp)).bad_ensembles);
                cur_sig=sig_pos(:,~obj.ensemble_filter(idx_ef(crp)).bad_ensembles);
                cur_vel=vel_data(:,~obj.ensemble_filter(idx_ef(crp)).bad_ensembles,:);
                cur_n=reshape(cur_n,[],1);
                cur_sig=reshape(cur_sig,[],1);
                cur_vel=reshape(cur_vel,[],3);
                fgood=isfinite(cur_n) & isfinite(cur_sig) & all(isfinite(cur_vel),2);
                cur_vel=cur_vel(fgood,:); % is this correct?
                idx=repmat(obj.mesh(idx_mesh(crp)).index(cur_n(fgood),cur_sig(fgood)),1,3);
                fgood=isfinite(idx);
                fcomp=cumsum(ones(size(idx)),2); % make an index to map component of velocity to output matrix
                vel{crp}=accumarray({idx(fgood),fcomp(fgood)},cur_vel(fgood),[obj.mesh(idx_mesh(crp)).ncells,3],@(x) mean(x,'all','omitnan'),nan,false); % compute mean velocity for each mesh cell
                velstd{crp}=accumarray({idx(fgood),fcomp(fgood)},cur_vel(fgood),[obj.mesh(idx_mesh(crp)).ncells,3],@(x) std(x,0,'all','omitnan'),nan,false); % compute standard deviation of velocity in each mesh cell
            end
        end
    end
end