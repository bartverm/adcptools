classdef CoherenceRegularization < Regularization
    properties
        C3
        C4
    end
    methods(Access = protected)
        function assemble_matrix_private(obj, solver)
            % C3 is the only matrix that can always be assembled
            % (independent of model)
            obj.C3 = assemble_coherence(obj, solver);

            assert(isa(solver.data_model,"TaylorModel"), "To constrain cell-based gradients, a Taylor model is required");
            if (all(solver.data_model.n_order > 0)) && (all(solver.data_model.sigma_order > 0) || all(solver.data_model.z_order > 0))
                % This condition may be relaxed a bit.
                obj.C4 = assemble_consistency(obj, solver);
            else
                warning('No consistency matrix assembled: Fit a sufficient number of Taylor terms')
            end
        end

        function C = assemble_coherence(obj, solver)

            Np = sum(solver.data_model.npars);
            Diag = speye(Np*solver.mesh.ncells);
            rows = []; cols = []; vals = [];
            for idx = 1:solver.mesh.ncells %rows
                nb = obj.neighbors(:, idx);
                nb = nb(~isnan(nb)); nnb = length(nb);
                for nb_idx = 1:nnb
                    col = ((nb(nb_idx)-1)*Np+1):(nb(nb_idx)*Np);
                    row = (idx-1)*Np+1:idx*Np;
                    val = -1/nnb*ones(1,Np);
                    rows = [rows row];
                    cols = [cols col];
                    vals = [vals val];
                end
            end
            C = Diag + sparse(rows, cols, vals, NNp, NNp);

            W = obj.assemble_weights(solver);

            C = W*C;

        end

        function C4 = assemble_consistency(obj, solver)

            % Matrix that matches cell-based gradients to inter-cell
            % gradients

            mesh = solver.mesh;
            wl = solver.water_level;

            rows = []; cols = []; vals = [];
            row_idx = 1;
            for idx = 1:mesh.ncells %rows
                for tid_idx = 1:length(wl.names)
                    tid_name = wl.names{tid_idx}(end-2:end); % Last three letters
                    [row, col, val] = obj.idx2element_mat(idx, tid_name, solver);
                    rows = [rows row_idx+row];
                    cols = [cols col];
                    vals = [vals val];
                    row_idx = max(rows) + 1;
                end
            end
            C4 = sparse(rows, cols, vals, 6*mesh.ncells*(2*length(solver.data_model.constituentsU)+1), mesh.ncells*sum(solver.data_model.npars));
        end

        function IM = assemble_incidence(obj, solver)
            mesh = solver.mesh;
            IM = zeros(mesh.ncells);
            for j = 1:mesh.ncells % loop trough every cell
                nbreal = obj.neighbors(j,~isnan(obj.neighbors(j,:)));
                nbreal = nbreal(nbreal>j);
                IM(j, nbreal) = 1;
            end
            IM = IM + IM';
        end

        function W = assemble_weights(obj, solver)
            par_names = solver.data_model.names;
            Np = sum(solver.data_model.npars);
            mesh = solver.mesh;
            w = ones([Np,1]);

            % Apply enhanced regularization for small singular value features using
            % characteristic spatial scales

            vertscale = max(mesh.z_patch) - min(mesh.z_patch);
            horscale = max(mesh.n_patch) - min(mesh.n_patch); %Typical scales

            % Automatically assign weights to smoothness of different parameters
            % Important due to orientation of main flow

            for i = 1:Np
                if contains(par_names{i}, 'w')
                    w(i) = w(i)*horscale/vertscale;
                end
                if contains(par_names{i}, 'dy') || contains(par_names{i}, 'dx')
                    w(i) = w(i)*horscale;
                end
            end
            wj = cell(1,mesh.ncells);
            wj(:)= {w};
            W = helpers.spblkdiag(wj{:});
        end

        function [row, col, val] = idx2element_mat(obj, idx, tid_name, solver)
            % Helper function for the extended consistency matrix assembly
            % (assembleC4.m)

            sig_center = solver.mesh.sig_center;
            n_center = solver.mesh.n_middle(mesh.col_to_cell);

            nb = obj.neighbors(:,idx);
            dom = obj.domains(idx);

            if dom==0 || dom==1 || dom == 5 % Internal cells (in terms of sigma)
                dsig = sig_center(nb(2))-sig_center(nb(4)); % central difference
            elseif dom ==2 || dom==3|| dom==4                           %Surface cells (deprecated: replaced by boundary condition matrix)
                dsig = sig_center(idx)-sig_center(nb(4)); % one-sided difference
            else                                                        % Bottom cells (deprecated: replaced by boundary condition matrix)
                dsig = sig_center(nb(2))-sig_center(idx); % one-sided difference
            end

            if dom==0 || dom==3 || dom == 7 % Internal cells (laterally)
                dn = n_center(nb(1))-n_center(nb(3)); % central difference
            elseif dom ==1 || dom==2|| dom==8                           % Right side of domain (deprecated)
                dn = n_center(idx)-n_center(nb(3)); % one-sided difference
            else                                                        % Left side of domain (deprecated)
                dn = n_center(nb(1))-n_center(idx); % one-sided difference
            end

            % First y derivatives, then sigma derivatives in order of u, v, w
            row = [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5]; % Row incremental index, always the same
            val = [1, -1/(dn), 1/(dn),...
                1, -1/(dn), 1/(dn),...
                1, -1/(dn), 1/(dn),...
                1, -1/(dsig), 1/(dsig),...
                1, -1/(dsig), 1/(dsig),...
                1, -1/(dsig), 1/(dsig)]; % Value, always the same (but look at order of magnitudes difference between the values...)

            pnames = obj.par_names_tot;
            col = [find(strcmp(pnames, sprintf('%s%i%s%s','cell ',idx,': d^1u/dy^1: ', tid_name))) , ...
                find(strcmp(pnames, sprintf('%s%i%s%s','cell ',nb(1),': u0: ', tid_name))) , ...
                find(strcmp(pnames, sprintf('%s%i%s%s','cell ',nb(3),': u0: ', tid_name))) , ...
                find(strcmp(pnames, sprintf('%s%i%s%s','cell ',idx,': d^1v/dy^1: ', tid_name))) , ...
                find(strcmp(pnames, sprintf('%s%i%s%s','cell ',nb(1),': v0: ', tid_name))) , ...
                find(strcmp(pnames, sprintf('%s%i%s%s','cell ',nb(3),': v0: ', tid_name))) , ...
                find(strcmp(pnames, sprintf('%s%i%s%s','cell ',idx,': d^1w/dy^1: ', tid_name))) , ...
                find(strcmp(pnames, sprintf('%s%i%s%s','cell ',nb(1),': w0: ', tid_name))) , ...
                find(strcmp(pnames, sprintf('%s%i%s%s','cell ',nb(3),': w0: ', tid_name))) , ...
                find(strcmp(pnames, sprintf('%s%i%s%s','cell ',idx,': d^1u/dsig^1: ', tid_name))) , ...
                find(strcmp(pnames, sprintf('%s%i%s%s','cell ',nb(2),': u0: ', tid_name))) , ...
                find(strcmp(pnames, sprintf('%s%i%s%s','cell ',nb(4),': u0: ', tid_name))) , ...
                find(strcmp(pnames, sprintf('%s%i%s%s','cell ',idx,': d^1v/dsig^1: ', tid_name))) , ...
                find(strcmp(pnames, sprintf('%s%i%s%s','cell ',nb(2),': v0: ', tid_name))) , ...
                find(strcmp(pnames, sprintf('%s%i%s%s','cell ',nb(4),': v0: ', tid_name))) , ...
                find(strcmp(pnames, sprintf('%s%i%s%s','cell ',idx,': d^1w/dsig^1: ', tid_name))) , ...
                find(strcmp(pnames, sprintf('%s%i%s%s','cell ',nb(2),': w0: ', tid_name))) , ...
                find(strcmp(pnames, sprintf('%s%i%s%s','cell ',nb(4),': w0: ', tid_name)))];

            % 18 terms in the sparse constructor element matrices. Reduce the number of
            % terms when not in the interior of the cross section.


            if dom == 1
                keep_idx = 10:18;
            elseif dom == 2
                keep_idx = [];
            elseif dom == 3
                keep_idx = 1:9;
            elseif dom == 4
                keep_idx = [];
            elseif dom == 5
                keep_idx = 10:18;
            elseif dom == 6
                keep_idx = [];
            elseif dom == 7
                keep_idx = 1:9;
            elseif dom == 8
                keep_idx = [];
            end
            col = col(keep_idx);
            row = row(keep_idx);
            val = val(keep_idx);
        end

    end
end



