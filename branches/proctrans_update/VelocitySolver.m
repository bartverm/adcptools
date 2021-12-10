classdef VelocitySolver < handle
    % Abstract base class to solve ADCP repeat transect velocity data
    %
    %   Subclasses should implement the get_solver_input method.
    %
    %   obj=VelocitySolver(...) allows to set properties upon construction.
    %   Depending on the class objects are assigned to a specific property:
    %   VMADCP -> obj.adcp (required)
    %   Mesh -> obj.mesh (required)
    %   Bathymetry -> obj.bathymetry
    %   XSection -> obj.xs
    %   Filter -> obj.ensemble_filter
    %   VelocityModel -> obj.velocity_model 
    %
    %  VelocitySolver properties:
    %   adcp - VMADCP object with the adcp data
    %   mesh - Mesh object defining a mesh
    %   bathy - Defines the bathymetry
    %   xs - Defines the cross-section
    %   ensemble_filter - defines repeat transects to be processed
    %   velocity_model - defines the velocity model to fit
    %
    %   VelocitySolver methods:
    %   get_parameters - get velocity model parameters on mesh
    %   get_velocity - get velocity on mesh
    %   rotate_to_xs - rotates velocity to cross-section direction
    %
    %   see also: VMADCP, Mesh, Bathymetry, XSection, Filter,
    %   TimeBasedVelocitySolver, LocationBasedVelocitySolver

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

        % VelocitySolver/velocity_model
        %
        %   Velocity model to use when solving for velocity
        %
        %   see also: VelocitySolver, VelocityModel
        velocity_model (1,1) VelocityModel
    end
    methods
        function obj=VelocitySolver(varargin)
            has_vmadcp=false;
            has_mesh=false;
            has_bathy=false;
            has_xs=false;
            has_model=false;
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
                elseif isa(cur_arg,'VelocityModel')
                    has_model=true;
                    obj.velocity_model=cur_arg;
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
            if ~has_model
                obj.velocity_model=VelocityModel;
            end
        end
    end
    methods (Abstract, Access=protected)
% Get input data for velocity solver
%       [vpos, vdat, xform] = get_solver_input(obj) returns the velocity
%       position, the velocity data, and the transformation matrix to get
%       from the velcity data to veloicity components in earth coordinates.
%
%       Subclasses should implement this function   
        [vpos, vdat, xform] = get_solver_input(obj)
    end
    methods

        function [vel, cov_vel, n_bvels] = get_velocity(obj)
            %   Get velocity from model parameters
            %
            %   vel=get_velocity(obj) computes the velocity in the mesh cells by using
            %   the velocity computed by combining the beam velocities measured at the
            %   same time (this is done by the ADCP class). Velocities that where
            %   measured within a mesh cell (this is determined by the Mesh class) are
            %   averaged.
            %
            %   [vel,std] = get_velocity(obj) also returns the standard
            %   deviation in the velocity
            [pars, cov_pars, n_bvels] = get_parameters(obj);
            [vel, cov_vel] = deal(cell(numel(pars,1)));
            idx_mesh = obj.make_indices();
            for crp = 1 : numel(pars)
                cmesh = obj.mesh(idx_mesh(crp));
                [vel{crp}, cov_vel{crp}] = ...
                    obj.velocity_model.get_velocity( ...
                    pars{crp}, ...
                    cov_pars{crp}, ...
                    cmesh.time);
            end
        end

        function [pars, cov_pars, n_vels]=get_parameters(obj)
            % Solve velocity model parameters
            %   
            %   [pars, cov_pars, n_vels]=get_parameters(obj) Obtain the
            %   parameters of the model by fitting it to the velocity data.
            %   cov_pars returns the covatiance matrix of the model
            %   parameters and n_vels returns the number of velocity
            %   components used for the computation.
            %
            % see also: VelocitySolver, get_velocity, Mesh, VMADCP

            % get velocity position, velocity data, and transformation
            % matrices to obtain earth velocity
            [vpos, vel_data, xform] = get_solver_input(obj);

            % indices to match repeat transects with corresponsing mesh
            [idx_mesh, idx_ef] = obj.make_indices();

            %%% compute s,n and sigma coordinates
            % get velocity position transverse to cross section
            [s_pos, n_pos] = obj.xs.xy2sn(vpos(:,:,:,1), vpos(:,:,:,2)); 

            % get bed elevation at velocity position. This could optionally
            % be done with the bed detection at the beam.
            zb_pos = obj.bathy.get_bed_elev(vpos(:,:,:,1), vpos(:,:,:,2)); 
            
            % compute sigma coordinate of velocities
            sig_pos = (vpos(:,:,:,3) - zb_pos) ./...
                (obj.adcp.water_level - zb_pos); 

            % Get time of velocities
            time = obj.adcp.time;

            % Initialize output
            pars = cell(numel(idx_ef),1);
            cov_pars = cell(numel(idx_ef),1);
            n_vels = cell(numel(idx_ef),1);

            % Loop over repeat transects
            for crp = 1:numel(idx_ef)

                % get selection of data of current repeat transect and
                % replicate singleton dimensions to get uniform dimensions
                % of inputs
                ens_filt = ~obj.ensemble_filter(idx_ef(crp)).bad_ensembles;
                cur_n = n_pos(:, ens_filt, :);
                cur_s = s_pos(:, ens_filt, :);
                cur_t = time(:, ens_filt,:);
                cur_t = repmat(cur_t, ...
                    size(vel_data, 1), 1, size(vel_data, 3));
                cur_sig = sig_pos(:, ens_filt, :);
                cur_vel = vel_data(:, ens_filt, :);
                cur_z = vpos(:, ens_filt, :, 3);
                cur_xform = xform(:, ens_filt, :, :);
                cur_xform = repmat(cur_xform,...
                    [size(vel_data, 1), 1, 1, 1]);

                % vectorize all input, at this point we don't care which
                % data were collected together
                cur_n = reshape(cur_n, [], 1); 
                cur_s = reshape(cur_s, [], 1);
                cur_t = reshape(cur_t, [], 1);
                cur_sig = reshape(cur_sig, [], 1);
                cur_vel = reshape(cur_vel, [], 1);
                cur_z = reshape(cur_z, [], 1);
                cur_xform = reshape(cur_xform, [], 3);

                % create filter selecting finite input
                fgood = find(isfinite(cur_n) & isfinite(cur_sig) &...
                    isfinite(cur_vel) & all(isfinite(cur_xform), 2));

                % find mesh cell input data belongs to
                cmesh = obj.mesh(idx_mesh(crp));
                cell_idx = cmesh.index(cur_n(fgood), cur_sig(fgood));

                % combine filter selecting finite input with filter
                % selecting data that belongs has a correspoding mesh cell
                % and filter input data accordingly
                fgood_idx = isfinite(cell_idx);
                cell_idx = cell_idx(fgood_idx);                
                fgood = fgood(fgood_idx);
                cur_vel = cur_vel(fgood);
                cur_xform = cur_xform(fgood,:);
                cur_n = cur_n(fgood);
                cur_s = cur_s(fgood);
                cur_z = cur_z(fgood);
                cur_t = cur_t(fgood);
                cur_sig = cur_sig(fgood);

                % compute velocity model input
                n_center = reshape(...
                    cmesh.n_middle(cmesh.col_to_cell), [], 1);
                cur_n = cur_n - n_center(cell_idx); % delta_n
                cur_z = cur_z - cmesh.z_center(cell_idx); % delta_z
                cur_sig = cur_sig - cmesh.sig_center(cell_idx); % delta_sig
                cur_t = cur_t - cmesh.time; % delta time

                % get model matrices
                [Mu, Mv, Mw] = obj.velocity_model.get_model(...
                    cur_t, cur_s, cur_n, cur_z, cur_sig);

                % combine velocity input to earth matrix with velocity 
                % model
                cur_xform = [...
                    Mu.*cur_xform(:,1),...
                    Mv.*cur_xform(:,2),... 
                    Mw.*cur_xform(:,3)];

                %%% combine velocity input with transformation matrices and
                %%% prepare data to be combined in one cell element per
                %%% mesh cell
                cur_dat = [cur_vel cur_xform];
                % number of columns in model matrix + one for velocity 
                % component
                ndat=size(cur_dat,2); 
                % make an index to map each model matrix parameter to a 
                % column in output cell array
                dat_idx = cumsum(ones(size(cell_idx, 1), ndat), 2); 
                cell_idx = repmat(cell_idx, ndat, 1); 
                cur_dat = reshape(cur_dat, [], 1); 
                dat_idx = reshape(dat_idx, [], 1);


                % gather data in cell array: every row corresponds to a
                % mesh cell, every column to a model parameter. First
                % column holds velocity component
                gather_dat = accumarray({cell_idx,dat_idx},... % indices
                    cur_dat,... % velocity and model matrices
                    [obj.mesh(idx_mesh(crp)).ncells,ndat],... % output size
                    @(x) {x},... % function collecting data in cell
                    {},... % value to fill when no data is avilable
                    false);  % output not sparse
                
                % reorganize gathered data into columns, which are
                % inputs to the model fitting function
                gather_dat = num2cell(gather_dat, 1);
                
                % fit model for each cell
                [t_pars, t_n_bvels, t_cov_pars] = cellfun( ...
                    @obj.fit_model, ... % function to fit model
                    gather_dat{:}, ... % input data for fitting
                    'UniformOutput', false); % output is non-scalar

                % replace empty cell elements with array with NaN values
                f_empty = cellfun(@isempty, t_pars);
                t_pars(f_empty) = {nan(ndat - 1, 1)};
                t_cov_pars(f_empty) = {nan(ndat - 1, ndat - 1)};
                

                %%% construct matrix from output cell array
                % reshape to allow subsequent concatenation
                t_cov_pars = cellfun( ...
                    @(x) shiftdim(x, -1), ...
                    t_cov_pars, ...
                    'UniformOutput',false);
                t_pars = cellfun( ...
                    @(x) shiftdim(x, -1), ...
                    t_pars, ...
                    'UniformOutput',false);
                pars{crp} = vertcat(t_pars{:});
                cov_pars{crp} = vertcat(t_cov_pars{:});
                n_vels{crp} = vertcat(t_n_bvels);
            end
        end



        function [vel, cov]=rotate_to_xs(obj, orig_vel, orig_cov)
        % Rotates velocity to direction of cross-section
        %
        %   vel=obj.rotate_to_xs(orig_vel) Rotates the velocity orig_vel
        %   to the direction of the cross-section.
        %
        %   [vel, cov]=obj.rotate_to_xs(orig_vel, orig_cov) also rotates
        %   the covariance matrix
        %
        % See also: VelocitySolver, get_velocity
            vel = cell(size(orig_vel));
            cov = cell(size(orig_cov));
            for crp = 1:numel(orig_vel)
                [vels, veln] = obj.xs.xy2sn_vel( ...
                    orig_vel{crp}(:, 1), ...
                    orig_vel{crp}(:, 2));
                vel{crp} = [vels, veln, orig_vel{crp}(:, 3)];
                if nargin > 2
                    Tsn = obj.xs.xy2sn_tens(orig_cov{crp}(:, 1:2, 1:2));
                    cov{crp} = cat(3, [Tsn, orig_cov{crp}(:,3,1:2)],...
                        orig_cov{crp}(:,:,3));
                end
            end
        end
    end
    methods(Static, Access=protected)
        function [model_pars, n_dat, cov_matrix] = fit_model(vel, varargin)
        % This function fits the given model to the existing data

            n_dat = numel(vel); % number of velocity data
            model_mat = [varargin{:}];  % model matrix
            n_pars = size(model_mat, 2); % number of model parameters
            % handle rank deficient model matrices
            if rank(model_mat) < n_pars
                model_pars = nan(n_pars, 1);
                cov_matrix = nan(n_pars, n_pars);
                return
            end
            % fit model tot data
            [model_pars, ~, ~, cov_matrix] = lscov(model_mat, vel);
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