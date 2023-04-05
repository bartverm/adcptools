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

        reg (1,1) Regularization

        opts (1,1) SolverOptions = SolverOptions
    end
    methods
        function obj=RegularizedSolver(varargin)
            obj = obj@Solver(varargin{:})
            % All abstract Solver properties have been included
            has_vmadcp = false;
            has_bathy = false;
            has_xs = false;
            has_opts = false;
            has_reg = false;
            for cnt_arg=1:nargin
                cur_arg=varargin{cnt_arg};
                if isa(cur_arg,'VMADCP')
                    has_vmadcp=true;
                    var = 'adcp';
                elseif isa(cur_arg, 'Bathymetry')
                    has_bathy=true;
                    continue
                elseif isa(cur_arg, 'Filter')
                    var = 'ensemble_filter';
                elseif isa(cur_arg, 'XSection')
                    has_xs=true;
                    continue
                elseif isa(cur_arg,'SolverOptions')
                    has_opts = true;
                    var = 'opts';
                elseif isa(cur_arg,'Regularization')
                    var = 'reg';
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

        function MP = get_parameters(obj)
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
            ds = s_pos;
            dn = n_pos - n_center(cell_idx); % delta_n
            dz = z_pos - cmesh.z_center(cell_idx); % delta_z
            dsig = sig_pos - cmesh.sig_center(cell_idx); % delta_sig

            % Data matrix: All data
            M = obj.data_model.get_model(time, ds, dn, dz, dsig);
%             obj.data_model.rotation_matrix
%             obj.data_model.rotation_matrix = inv( obj.data_model.rotation_matrix);
%             obj.data_model.rotation_matrix
%             M = obj.data_model.rotate_matrix(M);

            disp('Assembled model matrices')

            xform = xform*obj.data_model.rotation_matrix;
            Mb0 = [M(:,:,1).*xform(:,1),...
                M(:,:,2).*xform(:,2),...
                M(:,:,3).*xform(:,3)]; %Model matrix times unit vectors q
            disp('Assembled parameter - data mapping')
            [M, b, ns] = obj.reorder_model_matrix(Mb0, dat, cell_idx);
            
            MP = ModelParameters(M = M, b = b, reg = obj.reg, opts = obj.opts);

            MP.p = obj.solve(M,b);
            MP.ns = ns;
            disp('Finished')
        end

        function p = solve(obj, M, b)
            Np = size(M,2);
            % Generate first guess for p <-> estimate parameters
            p = nan([Np,length(obj.opts.reg_pars)]);
            Mg = M'*M;
            for idx = 1:length(obj.opts.reg_pars)
                rp = obj.opts.reg_pars{idx};
                Cg = sparse(0);
                for reg_idx = 1:numel(obj.reg.Cg)
                    Cg = Cg + rp(reg_idx)*obj.reg.Cg{reg_idx};
                end
                A = Mg + Cg; % Data plus constraints
                rhs = M'*b + rp(end)*obj.reg.C{end}'*obj.reg.rhs;
                [p(:,idx), flag] = obj.solve_single(A, rhs);
                disp(['Obtained solution using lambda = [', num2str(rp), ']^T after ', num2str(flag), ' iterations.'])
            end
        end

        function [p, iter] = solve_single(obj, A, rhs)
            if obj.opts.set_diagcomp
                obj.opts.preconditioner_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
            end
            L = ichol(A, obj.opts.preconditioner_opts);
            [p, ~, ~, iter, resvec] = pcg(A, rhs, obj.opts.pcg_tol, obj.opts.pcg_iter, L, L');
            
        end



        function [M, b, ns] = reorder_model_matrix(obj, Mb0, dat, cell_idx)
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
%             [vdat, xform] = obj.filter_and_vectorize(vdat, xform);

            vdat = obj.adcp.water_velocity(CoordinateSystem.Beam); % get velocity data

            % get transformation matrix
            xform = obj.adcp.xform(CoordinateSystem.Beam, CoordinateSystem.Earth); % get Earth to Beam transformation matrix

            xform(:,:,:,4) = []; % remove Error velocity to beam transformation
%             xform = obj.rotate_xform(xform);
            % filter and vectorize
            [vdat, xform] = obj.filter_and_vectorize(vdat, xform);

            
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

        function xform = rotate_xform(obj, xform)
            mat = [obj.xs.direction_orthogonal(1), obj.xs.direction_orthogonal(2), 0;
                   obj.xs.direction(1), obj.xs.direction(2), 0;
                   0, 0, 1];
%             xform(:,:,:,4)=[]; % remove Error velocity to beam transformation
            for i = 1:size(xform,2)
                xform(1,i,:,:) = squeeze(xform(1,i,:,:))*mat'; % Very shady, but good results!
            end

        end
    end

end