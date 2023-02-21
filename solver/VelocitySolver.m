classdef VelocitySolver < ADCPDataSolver
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

   methods
        function varargout = get_velocity(obj, varargin)
            %   Get velocity from model parameters
            %
            %   vel=get_velocity(obj) computes the velocity in the mesh cells by using
            %   the velocity computed by combining the beam velocities measured at the
            %   same time (this is done by the ADCP class). Velocities that where
            %   measured within a mesh cell (this is determined by the Mesh class) are
            %   averaged.
            %
            %   [vel,cov_vel] = get_velocity(obj) also returns the standard
            %   deviation in the velocity
            varargout = cell(1,nargout);
            [varargout{:}] = obj.get_data(varargin{:});
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
                cur_t = time(:, ens_filt, :);
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
                cur_t = seconds(cur_t);



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

        function [pars, cov_pars, n_vels]=get_parameters_reg(obj, opts)
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
                cur_t = time(:, ens_filt, :);
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
                dn = cur_n - n_center(cell_idx); % delta_n
                dz = cur_z - cmesh.z_center(cell_idx); % delta_z
                dsig = cur_sig - cmesh.sig_center(cell_idx); % delta_sig
                dt = cur_t - cmesh.time; % delta time
                dt = seconds(dt);
                %                 nb = length(dt);

                %                 [dt, cur_s, dn, dz, dsig,...
                %                     dtt, cur_st, dnt, dzt, dsigt] = split_dataset(opts, dt, cur_s, dn, dz, dsig);

                obj.velocity_model.get_parameter_names();

                % Internal continuity matrix
                %                 for j = 1:obj.mesh.ncells
                C1 = assembleC1(obj);
                %                 end
                C1p = C1'*C1;

                % External continuity matrix
                C2 = assembleC2(obj);
                C2p = C2'*C2;

                %Coherence matrix
                C3 = assembleC3(obj);
                C3p = C3'*C3;

                %Consistency matrix
                C4 = assembleC4(obj);
                C4p = C4'*C4;

                % Kinematic boundary condition matrix
                [C5, bc] = assembleC5(obj);
                C5p = C5'*C5;


                % Data matrix: All data
                [Mu, Mv, Mw] = obj.velocity_model.get_model(...
                    dt, cur_s, dn, dz, dsig);
                Mb0 = [Mu.*cur_xform(:,1), Mv.*cur_xform(:,2), Mw.*cur_xform(:,3)]; %Model matrix times unit vectors q

                % If cross-validation is to be applied: Split data in two
                % sets. All other operations before are data-independent
                % and can thus be performed only once.
                Mj = cell([obj.mesh.ncells,1]); ns = zeros([obj.mesh.ncells,1]); bj = cell([obj.mesh.ncells,1]);
                for j = 1:obj.mesh.ncells
                    Mj{j} = Mb0(j==cell_idx,:);
                    ns(j) = sum(j == cell_idx);
                    bj{j} = cur_vel(j==cell_idx);
                end

                b = vertcat(bj{:});
                M = spblkdiag(Mj{:});
                pcg_opts = struct('michol','on','type','ict','droptol',1e-3,'diagcomp',0);

                % Generate first guess for p
                rpg0 = [opts.reg_pars0{:}];
                A = M'*M + rpg0(1)*C1p + rpg0(2)*C2p + rpg0(3)*C3p + rpg0(4)*C4p + rpg0(5)*C5p;
                alpha = max(sum(abs(A),2)./diag(A))-2;
                pcg_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;                    
                L = ichol(A, pcg_opts);
                p0 = pcg(A, M'*b + rpg0(5)*C5'*bc, 1e-9, size(A,2), L, L'); % global, iterative inversion
                np = length(p0);

                % Compute quantities of interest for comparison
                rpg1 = opts.reg_pars1;
                p1 = nan(np, size(rpg1,1));
                for n = 1:length(rpg1{1})
                    A = M'*M + rpg1{1}(n)*C1p + rpg1{2}(n)*C2p + rpg1{3}(n)*C3p + rpg1{4}(n)*C4p + rpg1{5}(n)*C5p;
                    alpha = max(sum(abs(A),2)./diag(A))-2;
                    pcg_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;                   
                    p1(:,n) = pcg(A, M'*b + rpg1{5}(n)*C5'*bc, 1e-9, size(A,2), L, L'); % global, iterative inversion
                end



                % First loop: cross-validation ensembles
                opts.cell_idx = cell_idx;

                % Regularization parameters
                regP = combine_regpars(opts);

                nepochs = opts.cv_iter;
                p = zeros([np,size(regP,1),nepochs]);
                niter = nepochs*size(regP,1);
                i = 0;
                if opts.use_p0
                    pguess = p0;
                else
                    pguess = zeros(size(p0));
                end
                pe = zeros([10, size(regP,1), nepochs]);
                if ~strcmp(opts.cv_mode, 'none')
                    for ep = 1:nepochs
                        train_idx = logical(split_dataset(opts));
                        test_idx = ~train_idx;

                        % Construct training matrix and data
                        M0 = M(train_idx, :);
                        b0 = b(train_idx);

                        M1 = M(test_idx,:);
                        b1 = b(test_idx);

                        Mp = M0'*M0;

                        for rp = 1:size(regP,1)
                            i = i+1;
                            fprintf('Cross-validation percentage: %2.2f percent \n', 100*i/niter)
                            A = Mp + regP(rp,1)*C1p + regP(rp,2)*C2p + regP(rp,3)*C3p + regP(rp,4)*C4p + regP(rp,5)*C5p;
                            pcg_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
                            L = ichol(A, pcg_opts);
                            [p(:, rp, ep), ~, ~, it] = pcg(A, M0'*b0 + regP(rp,5)*C5'*bc, 1e-9, size(A,2), L, L', pguess); % Matrix of solutions (columns) belonging to regularization parameters regP (rows)

                            % Residuals and goodness of fit

                            pe(1, rp, ep) = calc_res(b, M*p(:, rp, ep)); % Performance on full set
                            pe(2, rp, ep) = calc_res(b0, M0*p(:, rp, ep)); % Performance on training set
                            pe(3, rp, ep) = calc_res(b1, M1*p(:, rp, ep)); % Performance on validation set
                            pe(4, rp, ep) = calc_res(0, C1*p(:, rp, ep)); % Performance on continuity
                            pe(5, rp, ep) = calc_res(0, C2*p(:, rp, ep)); % Performance on gen. continuity
                            pe(6, rp, ep) = calc_res(0, C3*p(:, rp, ep)); % Performance on smoothness
                            pe(7, rp, ep) = calc_res(0, C4*p(:, rp, ep)); % Performance on consistency
                            pe(8, rp, ep) = calc_res(bc, C5*p(:, rp, ep)); % Performance on boundary conditions

                            pe(9,rp,ep) = condest(A);
                            pe(10,rp,ep) = it;
                        end
                    end
                end
                Pe = mean(pe, 3);
                assignin("base", "dat", struct('M', M, 'C1', C1, 'C2', C2, 'C3', C3, 'C4', C4, 'C5', C5, 'bc', bc, ...
                    'IM', IM, 'p', p, 'p0', p0, 'p1', p1, 'b', b,...
                    'opts', opts, 'cell_idx', cell_idx, 'Pe', Pe, 'regP', regP))

                pars{1,1} = reshape(squeeze(p(:,1,1)) ,[size(Mb0,2), obj.mesh.ncells])'; % At this stage, pars are already known.
                cov_pars{1,1} = 0; n_vels{1,1} = ns;
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
        % Rotates velocity to direction of cross-section
        %
        %   vel=obj.rotate_to_xs(orig_vel) Rotates the velocity orig_vel
        %   to the direction of the cross-section.
        %
        %   [vel, cov]=obj.rotate_to_xs(orig_vel, orig_cov) also rotates
        %   the covariance matrix
        %
        % See also: VelocitySolver, get_velocity
            xs = [obj.xs];
            u = cellfun(@(x) x(:,1),orig_vel,'UniformOutput',false);
            v = cellfun(@(x) x(:,2),orig_vel,'UniformOutput',false);
            [vels, veln] = xs.xy2sn_vel( ...
                u, ...
                v);
            vel = cellfun(@(a,b,c) [a,b,c(:,3)], vels, veln, orig_vel, ...
                'UniformOutput',false);
            if nargin > 2
                Trot = cellfun(@(x) x(:,1:2,1:2), orig_cov, ...
                    'UniformOutput',false);
                Tsn = xs.xy2sn_tens(Trot);
                cov = cellfun(@(a,b) cat(3, [a, b(:,3,1:2)],...
                    b(:,:,3)),Tsn,orig_cov, ...
                    'UniformOutput',false);
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
            assignin("base","Mbj0",model_mat)
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