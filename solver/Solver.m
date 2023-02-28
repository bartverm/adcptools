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
        data_model (1,1) DataModel = DataModel;

        water_level (1,1) WaterLevel = WaterLevel;

        opts (1,1) SolverOptions = SolverOptions;

        regularizations (1,:) Regularization = [ContinuityRegularization; 
            CoherenceRegularization;
            KinematicRegularization]
    end
    methods
        function obj=Solver(varargin)
            obj = obj@helpers.ArraySupport(varargin{:})
            for cnt_arg=1:nargin
                cur_arg=varargin{cnt_arg};
                if isa(cur_arg, 'Mesh')
                    var = 'mesh';
                elseif isa(cur_arg, 'Bathymetry')
                    var = 'bathy';
                elseif isa(cur_arg,'XSection')
                    var = 'xs';
                elseif isa(cur_arg,'DataModel')
                    var = 'data_model';
                elseif isa(cur_arg,'SolverOptions')
                    var = 'opts';
                elseif isa(cur_arg,'WaterLevel')
                    var = 'water_level';
                else
                    continue
                end
                obj.assign_property(var,cur_arg);
            end
        end

        function [pars, cov_pars, n_vels] = get_parameters(obj)
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
                if strcmp(obj.opts.algorithm, 'lscov')
                    [pars, cov_pars, n_vels] = obj.run_method('get_parameters');
                    return
                elseif strcmp(obj.opts.algorithm, 'pcg')
                    [pars, cov_pars, n_vels] = obj.run_method('get_parameters_reg');
                    return
                end
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

        function [varargout] = get_data(obj, varargin)
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
            if nargin == 1
                varargin = cell(1,3);
                [varargin{:}] = get_parameters(obj);
            end
            if nargin > 1
                if nargin - 1 < nargout
                    error('VelocitySolver:WrongInputNumber',...
                        'Please pass no input other than object, or as many inputs as required outputs')
                end
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

        function mp = get_parameters_reg(obj)
            % Solve velocity model parameters
            %
            %   [pars, cov_pars, n_vels]=get_parameters_reg(obj) Obtain the
            %   parameters of the model by fitting it to the velocity data.
            %   cov_pars returns the covatiance matrix of the model
            %   parameters and n_vels returns the number of velocity
            %   components used for the computation.
            %
            % see also: Solver, get_velocity, Mesh, VMADCP

            if ~isscalar(obj)
                mp = obj.run_method('get_parameters_reg');
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
            ds = s_pos;
            dn = n_pos - n_center(cell_idx); % delta_n
            dz = z_pos - cmesh.z_center(cell_idx); % delta_z
            dsig = sig_pos - cmesh.sig_center(cell_idx); % delta_sig
            dt = time - cmesh.time; % delta time
            dt = seconds(dt);
            

            obj.velocity_model.get_parameter_names();

            all_C=obj.regularizations.assemble_matrix;

            % Internal continuity matrix
            %                 for j = 1:obj.mesh.ncells
            C1 = assembleC1(obj);
            %                 end
            C1p = C1'*C1;

            % External continuity matrix
            C2 = assembleC2(obj);
            C2p = C2'*C2;

            %Coherence matrix
            [C3, IM] = assembleC3(obj);
            C3p = C3'*C3;

            %Consistency matrix
            C4 = assembleC4(obj);
            C4p = C4'*C4;

            % Kinematic boundary condition matrix
            [C5, bc] = assembleC5(obj);
            C5p = C5'*C5;

            % Data matrix: All data
            [Mu, Mv, Mw] = obj.velocity_model.get_model(...
                dt, ds, dn, dz, dsig);
            Mb0 = [Mu.*xform(:,1), Mv.*xform(:,2), Mw.*xform(:,3)]; %Model matrix times unit vectors q

            Mj = cell([obj.mesh.ncells,1]); ns = zeros([obj.mesh.ncells,1]); bj = cell([obj.mesh.ncells,1]);
            for j = 1:obj.mesh.ncells
                Mj{j} = Mb0(j==cell_idx,:);
                ns(j) = sum(j == cell_idx);
                bj{j} = dat(j==cell_idx);
            end

            b = vertcat(bj{:});
            M = spblkdiag(Mj{:});

            mp.model = obj.velocity_model;
            mp.M = M;
            mp.C = {C1, C2, C3, C4, C5};
            mp.bc = bc;
            
            Np = size(M,2);
            % Generate first guess for p <-> estimate parameters
            p = nan([Np,length(obj.opts.reg_pars)]);
            for idx = 1:length(obj.opts.reg_pars)
                rp = obj.opts.reg_pars{idx};
                A = M'*M + rp(1)*C1p + rp(2)*C2p + rp(3)*C3p + rp(4)*C4p + rp(5)*C5p;
                if obj.opts.set_diagcomp
                    obj.opts.preconditioner_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
                end
                L = ichol(A, obj.opts.preconditioner_opts);
                p(:,idx) = pcg(A, M'*b + rp(5)*C5'*bc, obj.opts.pcg_tol, obj.opts.pcg_iter, L, L');
            end
            mp.p = p;
            

            % If cross-validation is to be applied: Split data in two
            % sets. All other operations before are data-independent
            % and can thus be performed only once.

            if ~strcmp(obj.cv_mode, 'random')
                % Only one partition involved in cv analysis
                train_idx = logical(obj.split_dataset(cell_idx));
            else
                train_idx = nan(length(cell_idx), obj.opts.cv_iter);
                for cvi = 1:obj.opts.cv_iter
                    train_idx(:,cvi) = logical(obj.split_dataset(cell_idx));
                end
                
            end
            test_idx = ~train_idx;


            for idx = 1:length(obj.opts.reg_pars)
                for tidx = 1:size(train_idx,2)
                    M_train = M(train_idx(:,tidx), :);
                    b_train = b(train_idx(:,tidx));

                    M_test = M(test_idx(:,tidx),:);
                    b_test = b(test_idx(:,tidx));

                    Mp_train = M_train'*M_train;
                    rp = obj.opts.reg_pars{idx};
                    A = Mp_train + rp(1)*C1p + rp(2)*C2p + rp(3)*C3p + rp(4)*C4p + rp(5)*C5p;
                    if obj.opts.set_diagcomp
                        obj.opts.preconditioner_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
                    end
                    L = ichol(A, obj.opts.preconditioner_opts);
                    p(:,idx) = pcg(A, M_train'*b_train + rp(5)*C5'*bc, obj.opts.pcg_tol, obj.opts.pcg_iter, L, L');
            
            end
% 
% 
% 
%             % First loop: cross-validation ensembles
%             opts.cell_idx = cell_idx;
% 
%             % Regularization parameters
%             
% 
%             nepochs = opts.cv_iter;
%             p = zeros([np,size(regP,1),nepochs]);
%             niter = nepochs*size(regP,1);
%             i = 0;
%             if opts.use_p0
%                 pguess = p0;
%             else
%                 pguess = zeros(size(p0));
%             end
% 
% 
% 
% 
% 
%             if obj.opts.gen_analysis
% 
%             pe = zeros([10, size(regP,1), nepochs]);
%             if ~strcmp(opts.cv_mode, 'none')
%                 for ep = 1:nepochs
%                     train_idx = logical(split_dataset(opts));
%                     test_idx = ~train_idx;
% regP = combine_regpars(opts);
%                     % Construct training matrix and data
%                     M0 = M(train_idx, :);
%                     b0 = b(train_idx);
% 
%                     M1 = M(test_idx,:);
%                     b1 = b(test_idx);
% 
%                     Mp = M0'*M0;
% 
%                     for rp = 1:size(regP,1)
%                         i = i+1;
%                         fprintf('Cross-validation percentage: %2.2f percent \n', 100*i/niter)
%                         A = Mp + regP(rp,1)*C1p + regP(rp,2)*C2p + regP(rp,3)*C3p + regP(rp,4)*C4p + regP(rp,5)*C5p;
%                         pcg_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
%                         L = ichol(A, pcg_opts);
%                         [p(:, rp, ep), ~, ~, it] = pcg(A, M0'*b0 + regP(rp,5)*C5'*bc, 1e-9, size(A,2), L, L', pguess); % Matrix of solutions (columns) belonging to regularization parameters regP (rows)
% 
%                         % Residuals and goodness of fit
% 
%                         pe(1, rp, ep) = calc_res(b, M*p(:, rp, ep)); % Performance on full set
%                         pe(2, rp, ep) = calc_res(b0, M0*p(:, rp, ep)); % Performance on training set
%                         pe(3, rp, ep) = calc_res(b1, M1*p(:, rp, ep)); % Performance on validation set
%                         pe(4, rp, ep) = calc_res(0, C1*p(:, rp, ep)); % Performance on continuity
%                         pe(5, rp, ep) = calc_res(0, C2*p(:, rp, ep)); % Performance on gen. continuity
%                         pe(6, rp, ep) = calc_res(0, C3*p(:, rp, ep)); % Performance on smoothness
%                         pe(7, rp, ep) = calc_res(0, C4*p(:, rp, ep)); % Performance on consistency
%                         pe(8, rp, ep) = calc_res(bc, C5*p(:, rp, ep)); % Performance on boundary conditions
% 
%                         pe(9,rp,ep) = condest(A);
%                         pe(10,rp,ep) = it;
%                     end
%                 end
%             end
%             Pe = mean(pe, 3);
%             assignin("base", "dat", struct('M', M, 'C1', C1, 'C2', C2, 'C3', C3, 'C4', C4, 'C5', C5, 'bc', bc, ...
%                 'IM', IM, 'p', p, 'p0', p0, 'p1', p1, 'b', b,...
%                 'opts', opts, 'cell_idx', cell_idx, 'Pe', Pe, 'regP', regP))
% 
%             pars{1,1} = reshape(squeeze(p(:,1,1)) ,[size(Mb0,2), obj.mesh.ncells])'; % At this stage, pars are already known.
%             cov_pars{1,1} = 0; n_vels{1,1} = ns;
         end
         end

        function training_idx = split_dataset(obj, cell_idx)

            training_idx = ones(size(cell_idx));

            if strcmp(obj.opts.cv_mode, 'none')
                training_idx = ones(size(cell_idx));
            elseif strcmp(obj.opts.cv_mode, 'random')
                tp = obj.opts.training_perc;
                rand0 = rand(size(training_idx));
                training_idx = (rand0 <= tp);
            elseif strcmp(obj.opts.cv_mode, 'omit_cells')
                oc = obj.opts.omit_cells;
                for occ = 1:length(oc)
                    training_idx(cell_idx==oc(occ)) = 0;
                end
            elseif strcmp(obj.opts.cv_mode, 'omit_time') % to be implemented
            end
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