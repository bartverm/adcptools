classdef Solver < helpers.ArraySupport
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
        xs (1,1) XSection = XSection;

        % Solver/data_model
        %
        %   Velocity model to use when solving for velocity
        %
        %   see also: Solver, DataModel
        data_model (1,1) DataModel = VelocityModel;

        % Solver/regularization
        %
        %   Regularization object. If not provided upon construction, an
        %   empty object will be created.
        %
        %   see also: Solver, DataModel
        regularization (1,:) regularization.Regularization =...
            regularization.Regularization


        % Solver/opts
        %
        %   SolverOptions object
        %
        %   see also: Solver, DataModel
        opts (1,1) SolverOptions = SolverOptions

    end
    methods
        % function obj=Solver(varargin)
        %     obj = obj@helpers.ArraySupport(varargin{:})
        %     obj.parse_class_params_inputs(varargin{:});
        % end

        function MP = get_parameters(obj)
            % Solve velocity model parameters
            %
            %   [pars, cov_pars, n_vels]=get_parameters(obj) Obtain the
            %   parameters of the model by fitting it to the velocity data.
            %   cov_pars returns the covatiance matrix of the model
            %   parameters and n_vels returns the number of velocity
            %   components used for the computation.
            %
            % see also: Solver, get_velocity, Mesh, VMADCP

            if ~isscalar(obj)
                MP = obj.run_method('get_parameters');
                return
            end

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
            d_s = s_pos;
            d_n = n_pos - n_center(cell_idx); % delta_n
            d_z = z_pos - cmesh.z_center(cell_idx); % delta_z
            d_sig = sig_pos - cmesh.sig_center(cell_idx); % delta_sig

            % get model matrices
            M = obj.data_model.get_model(...
                time, d_s, d_n, d_z, d_sig);

            disp('Assembled model matrices')
            % From here, different solvers are used depending on the
            % SolverOptions opts

            switch obj.opts.algorithm
                case "lscov"
                    % combine velocity input to earth matrix with velocity
                    % model
                    npars = obj.data_model.npars;
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
                    MP = ModelParameters(opts = obj.opts);
                    MP.pars = vertcat(t_pars{:});
                    MP.cov_pars = vertcat(t_cov_pars{:});
                    MP.ns = vertcat(t_n_bvels{:});
                    MP.p = MP.pars2p();
                case "pcg"
                    npars=obj.data_model.npars;
                    Mb0 = [M(:,1:npars(1),1).*xform(:,1),...
                        M(:,1:npars(2),2).*xform(:,2),...
                        M(:,1:npars(3),3).*xform(:,3)]; %Model matrix times unit vectors q
                    disp('Assembled parameter - data mapping')
                    [M, b, ns] = obj.reorder_model_matrix(Mb0, dat, cell_idx);
                    MP = ModelParameters(M = M, b = b, opts = obj.opts,...
                        regularization = obj.regularization);
                    MP.p = obj.solve(M,b);
                    MP.ns = ns;
                    disp('Finished')
                    MP.pars = MP.p2pars();
                otherwise
                    error("Enter correct algorithm in SolverOptions: lscov (only data, cell-based) or pcg (regularized, global)")
            end
        end

        function [varargout] = get_data(obj, mp, varargin)
            %   Get data from model parameters
            %
            %   vel=get_data(obj)
            %
            %   [vel,cov_vel] = get_data(obj)

            % Handle non-scalar call
            varargout = cell(1,nargout);
            if ~isscalar(obj)
                [varargout{:}] = obj.run_method('get_data',varargin{:});
                return
            end

            % scalar call
            if nargin < 2
                mp = get_parameters(obj);
            end
            varargout = cell(1,nargout);
            [varargout{:}] = ...
                obj.data_model.get_data(varargin{:});
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

        function p = solve(obj, M, b)
            Np = size(M,2);
            % Generate first guess for p <-> estimate parameters
            reg_pars = {obj.regularization.weight};
            siz_reg = cellfun(@size,reg_pars,'UniformOutput',false);
            assert(isscalar(reg_pars) || isequal(siz_reg{:}),...
                'Weights of regression objects must have the same size')
            reg_pars = vertcat(reg_pars{:});
            n_sols = size(reg_pars,2);
            n_regs = size(reg_pars,1);
            p = nan([Np,n_sols]);
            Mg = M'*M; % Expensive operation -> minimize number of calls
            for idx = n_sols
                rp = reg_pars(:,idx);
                Cg = sparse(0);
                rhs = M'*b;
                for reg_idx = 1:n_regs
                    Cg = Cg + rp(reg_idx)*obj.regularization(reg_idx).Cg;
                    rhs = rhs + ...
                        rp(reg_idx)*obj.regularization(reg_idx).C'*...
                        obj.regularization(reg_idx).rhs;
                end
                A = Mg + Cg; % Data plus constraints
                [p(:,idx), flag] = obj.solve_single(A, rhs);
                disp(['Obtained solution using lambda = [', num2str(rp'), ']^T after ', num2str(flag), ' iterations.'])
            end
        end

        function [p, iter] = solve_single(obj, A, rhs)
            if obj.opts.set_diagcomp
                obj.opts.preconditioner_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
            end
            L = ichol(A, obj.opts.preconditioner_opts);
            [p, ~, ~, iter, ~] = pcg(A, rhs, obj.opts.pcg_tol, obj.opts.pcg_iter, L, L');

        end



        function [M, b, ns] = reorder_model_matrix(obj, Mb0, dat, cell_idx)
            ncells = obj.mesh.ncells;
            npars = size(Mb0,2);
            nbvels = size(Mb0,1);

            % sort cells, to diaganolize matrix
            [cell_idx, srt_idx] = sort(cell_idx);
            Mb0 = Mb0(srt_idx,:);
            b = sparse(dat(srt_idx));

            % build sparse matrix indices
            col_idx = (cell_idx - 1) * npars + (1 : npars);
            row_idx = repmat((1 : nbvels)', 1, npars);
            M = sparse(row_idx(:), col_idx(:), Mb0(:),...
                nbvels, npars * ncells);
            ns = accumarray(cell_idx,ones(size(cell_idx)),[ncells, 1]);
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
end