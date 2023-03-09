classdef RegularizedSolver < Solver
    properties
        % Solver/adcp
        %
        %   Scalar VMADCP object holding the adcp data to compute the velocity
        %
        %   see also: Solver, VMADCP
        adcp (:,1) VMADCP {mustBeScalarOrEmpty} = rdi.VMADCP.empty

        % Solver/ensemble_filter
        %
        %   Ensemble filter defining ensembles to include
        %
        %   see also: Solver, EnsembleFilter
        ensemble_filter (1,1) EnsembleFilter

        regularization (1,1) Regularization

        opts (1,1) SolverOptions
    end
    methods
        function obj=RegularizedSolver(varargin)
            obj = obj@Solver(varargin{:})
            has_vmadcp=false;
            has_bathy=false;
            has_xs=false;
            for cnt_arg=1:nargin
                cur_arg=varargin{cnt_arg};
                if isa(cur_arg,'VMADCP')
                    has_vmadcp=true;
                    var = 'adcp';
                elseif isa(cur_arg, 'Bathymetry')
                    has_bathy=true;
                    continue
                elseif isa(cur_arg,'Filter')
                    var = 'ensemble_filter';
                elseif isa(cur_arg,'XSection')
                    has_xs=true;
                elseif isa(cur_arg,'SolverOptions')
                    var = 'opts';
                elseif isa(cur_arg,'Regularization')
                    var = 'regularization';
                    has_reg=true;
                else
                    continue
                end
                obj.assign_property(var, cur_arg)
            end

            if ~has_bathy && has_vmadcp
                B=BathymetryScatteredPoints(obj.adcp);
                obj.assign_property('bathy',B);
            end
            if ~has_xs && has_vmadcp
                XS=XSection(obj.adcp);
                obj.assign_property('xs',XS);
            end

            if ~has_reg % Order of inserting properties is important
                R = Regularization(bathy = obj.bathy, xs = obj.xs, mesh = obj.mesh, model = obj.data_model);
                R.assemble_matrices()
                obj.assign_property('regularization', R);
            end
        end

        function mp = get_parameters(obj)
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
                mp = obj.run_method('get_parameters');
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

            % Data matrix: All data
            M = obj.data_model.get_model(...
                dt, ds, dn, dz, dsig);

            Mb0 = [M(:,:,1).*xform(:,1),...
                M(:,:,2).*xform(:,2),...
                M(:,:,3).*xform(:,3)]; %Model matrix times unit vectors q

            [M, b] = obj.reorder_model_matrix(Mb0, dat, cell_idx);
            MP = ModelParameters(M = M, b = b, reg = obj.regularization, opts = obj.opts);
            mp.model = obj.velocity_model;
            mp.M = M;
            mp.C = {C1, C2, C3, C4, C5};
            mp.bc = bc;

            Np = size(M,2);
            % Generate first guess for p <-> estimate parameters
            p = nan([Np,length(obj.opts.reg_pars)]);
            Mp.p = obj.solve()
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
        end

        function p = solve()
            % Continue here!
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

        function [M, b] = reorder_model_matrix(obj, Mb0, dat, cell_idx)
            Mj = cell([obj.mesh.ncells,1]); 
            ns = zeros([obj.mesh.ncells,1]); 
            bj = cell([obj.mesh.ncells,1]);
            for j = 1:obj.mesh.ncells
                Mj{j} = Mb0(j==cell_idx,:);
                ns(j) = sum(j == cell_idx);
                bj{j} = dat(j==cell_idx);
            end
            b = vertcat(bj{:});
            M = helpers.spblkdiag(Mj{:});
        end

    end



    methods (Access=protected)
        function [vpos, vdat, xform, time, wl] = get_solver_input(obj)
            vpos = obj.adcp.depth_cell_position; % velocity positions
            vdat = [];
            xform = [];
            time = obj.adcp.time;
            wl = obj.adcp.water_level;

            % get selection of data of current repeat transect and
            % replicate singleton dimensions to get uniform dimensions
            % of inputs
            ens_filt = ~obj.ensemble_filter.bad_ensembles;
            vpos = vpos(:, ens_filt, :,:);
            time = time(ens_filt);
            time = repmat(time, size(vpos, 1), 1, size(vpos, 3));
            wl = wl(ens_filt);
            wl = repmat(wl, size(vpos, 1), 1, size(vpos, 3));

            % vectorize
            time = reshape(time, [], 1);
            vpos = reshape(vpos, [], 3);
            wl = reshape(wl, [], 1);
        end

        function [vdat,xform] = filter_and_vectorize(obj,vdat, xform)
            ens_filt = ~obj.ensemble_filter.bad_ensembles;

            % ensemble filter
            xform = xform(:, ens_filt, :, :);
            xform = repmat(xform,...
                [size(vdat, 1), 1, 1, 1]);
            vdat = vdat(:,ens_filt,:);

            % vectorize
            vdat = reshape(vdat, [], 1);
            xform = reshape(xform, [], size(xform,4));
        end
    end

end