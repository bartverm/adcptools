classdef Solver < handle
    % Abstract base class to solve ADCP repeat transect velocity data
    %
    %   Subclasses should implement the get_solver_input method.
    %
    %   obj=Solver(...) allows to set properties upon construction.
    %   Depending on the class objects are assigned to a specific property:
    %   VMADCP -> obj.adcp (required)
    %   Mesh -> obj.mesh (required)
    %   Bathymetry -> obj.bathymetry
    %   XSection -> obj.xs
    %   Filter -> obj.ensemble_filter
    %   DataModel -> obj.data_model
    %
    %  Solver properties:
    %   adcp - VMADCP object with the adcp data
    %   mesh - Mesh object defining a mesh
    %   bathy - Defines the bathymetry
    %   xs - Defines the cross-section
    %   ensemble_filter - defines repeat transects to be processed
    %   data_model - defines the velocity model to fit
    %
    %   Solver methods:
    %   get_parameters - get velocity model parameters on mesh
    %   get_velocity - get velocity on mesh
    %   rotate_to_xs - rotates velocity to cross-section direction
    %
    %   see also: VMADCP, Mesh, Bathymetry, XSection, Filter,
    %   VelocitySolver, BackscatterSolver


    properties
        % Solver/mesh mesh on which velocity is solved
        %
        %   Mesh object defining the mesh on which the data is to be
        %   solved. If an array is passed, each mesh is used for a
        %   different repeat transect, thus the number of elements must
        %   match the number of elements in the ensemble_filter property.
        %
        %   Default value is SigmaZetaMesh
        %
        %   see also: Solver, Mesh
        mesh (1,1) Mesh = SigmaZetaMesh;

        % Solver/bathy
        %
        %   Scalar bathymetry object that defines the location of the bed. Default
        %   is BathymetryFromScatteredPoints(adcp), which constructs a bathymetry
        %   based on the VMADCP data in adcp.
        %
        %   see also: Solver, BathymetryScatteredPoints
        bathy (1,1) Bathymetry = BathymetryScatteredPoints;


        % Solver/xs
        %
        %   Scalar XSection object defining the origin and direction of the
        %   cross-section. Default value is XSection(adcp) which construct a
        %   cross-section based on the track of the VMADCP data in adcp
        %
        %   see also: Solver, XSection
        xs (1,1) XSection

        % Solver/data_model
        %
        %   Velocity model to use when solving for velocity
        %
        %   see also: Solver, DataModel
        data_model (1,1) DataModel
    end
    methods
        function obj=Solver(varargin)
            has_mesh=false;
            has_model=false;
            for cnt_arg=1:nargin
                cur_arg=varargin{cnt_arg};
                if isa(cur_arg, 'Mesh')
                    has_mesh=true;
                    obj.mesh=cur_arg;
                elseif isa(cur_arg, 'Bathymetry')
                    obj.bathy=cur_arg;
                elseif isa(cur_arg,'XSection')
                    obj.xs=cur_arg;
                elseif isa(cur_arg,'DataModel')
                    has_model=true;
                    obj.data_model=cur_arg;
                end
            end
            if ~has_mesh
                error('You must provide a Mesh object upon construction of a Solver object')
            end
            if ~has_model
                obj.data_model=DataModel;
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
            % see also: Solver, get_velocity, Mesh, VMADCP

            % get velocity position, velocity data, and transformation
            % matrices to obtain earth velocity
            [vpos, dat, xform, time, wl] = get_solver_input(obj);

            %%% compute s,n and sigma coordinates
            % get velocity position transverse to cross section
            [s_pos, n_pos] = obj.xs.xy2sn(vpos(:,1), vpos(:,2));

            % get bed elevation at velocity position. This could optionally
            % be done with the bed detection at the beam.
            zb_pos = obj.bathy.get_bed_elev(vpos(:,1), vpos(:,2));

            % compute sigma coordinate of velocities
            z_pos=vpos(:,3);
            sig_pos = (z_pos - zb_pos) ./...
                (wl - zb_pos);

            % create filter selecting finite input
            fgood = find(isfinite(n_pos) & isfinite(sig_pos) &...
                isfinite(dat) & all(isfinite(xform), 2));

            % find mesh cell input data belongs to
            cmesh = obj.mesh;
            cell_idx = cmesh.index(n_pos(fgood), sig_pos(fgood));

            % combine filter selecting finite input with filter
            % selecting data that belongs has a correspoding mesh cell
            % and filter input data accordingly
            fgood_idx = isfinite(cell_idx);
            cell_idx = cell_idx(fgood_idx);
            fgood = fgood(fgood_idx);
            dat = dat(fgood);
            xform = xform(fgood,:);
            n_pos = n_pos(fgood);
            s_pos = s_pos(fgood);
            z_pos = z_pos(fgood);
            time = time(fgood);
            sig_pos = sig_pos(fgood);

            % compute velocity model input
            n_center = reshape(...
                cmesh.n_middle(cmesh.col_to_cell), [], 1);
            n_pos = n_pos - n_center(cell_idx); % delta_n
            z_pos = z_pos - cmesh.z_center(cell_idx); % delta_z
            sig_pos = sig_pos - cmesh.sig_center(cell_idx); % delta_sig
            time = time - cmesh.time; % delta time
            time = seconds(time);

            % get model matrices
            M = obj.data_model.get_model(...
                time, s_pos, n_pos, z_pos, sig_pos);
            npars = obj.data_model.npars;

            % combine velocity input to earth matrix with velocity
            % model
            assert(size(xform,2)<=size(M,3),...
                'Solver:WrongNComponents',...
                'More component in transformation matrix than in model matrix')
            
            ncomp = size(xform,2); % number of component in transformation matrix
            npars(ncomp+1:end)=[]; % only keep number of parameters for component in transformation matrix
            xformM = nan(size(xform,1), sum(npars));
            cum_pars= cumsum([0 npars]);
            for ccomp = 1 : size(xform,2)
                xformM(:,cum_pars(ccomp)+1:cum_pars(ccomp+1)) = ...
                    M(:,1:npars(ccomp),ccomp).*xform(:,ccomp);
            end
            xform = xformM;

            %%% combine velocity input with transformation matrices and
            %%% prepare data to be combined in one cell element per
            %%% mesh cell
            dat = [dat xform];
            % number of columns in model matrix + one for velocity
            % component
            ndat=size(dat,2);
            % make an index to map each model matrix parameter to a
            % column in output cell array
            dat_idx = cumsum(ones(size(cell_idx, 1), ndat), 2);
            cell_idx = repmat(cell_idx, ndat, 1);
            dat = reshape(dat, [], 1);
            dat_idx = reshape(dat_idx, [], 1);


            % gather data in cell array: every row corresponds to a
            % mesh cell, every column to a model parameter. First
            % column holds velocity component
            gather_dat = accumarray({cell_idx,dat_idx},... % indices
                dat,... % velocity and model matrices
                [obj.mesh.ncells,ndat],... % output size
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
            pars = vertcat(t_pars{:});
            cov_pars = vertcat(t_cov_pars{:});
            n_vels = vertcat(t_n_bvels{:});
        end
    end
    methods (Abstract, Access=protected)
        % Get input data for velocity solver
        %       [vpos, vdat, xform, time] = get_solver_input(obj) returns the velocity
        %       position, the velocity data, the transformation matrix  to get
        %       from the velocity data to velocity components in earth
        %       coordinates and the time of the data
        %
        %       Subclasses should implement this function
        [vpos, vdat, xform, time, wl] = get_solver_input(obj)
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