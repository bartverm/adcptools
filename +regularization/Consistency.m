classdef Consistency < regularization.TaylorBased
    methods(Access = protected)
        function assemble_matrix_private(obj)
            obj.mustBeTaylorModel;

            % find n and s coordinate of cells
            n_c = obj.mesh.n_middle(obj.mesh.col_to_cell);
            sig_c = reshape(obj.mesh.sig_center,1,[]);

            % make indices and values
            [row_n, col_n, val_n] = obj.make_idx_val('n', [1 3], n_c);
            [row_s, col_s, val_s] = obj.make_idx_val('sig', [2 4], sig_c);
            row_s = row_s + size(row_n,1);

            % build matrix
            npars = sum(obj.model.npars);
            ncells = obj.mesh.ncells;
            neq = size(row_n,1) + size(row_s,1);
            obj.C = sparse(...
                [row_n(:); row_s(:)],...
                [col_n(:); col_s(:)],...
                [val_n(:); val_s(:)],...
                neq,...
                ncells * npars);
        end

    end
    methods(Access = protected)
        function val = get_min_order(obj)
            nc = obj.model.ncomponents;
            val = zeros([5,nc]);
        end
        function [row, col, val] = make_idx_val(obj, var_name, nb_idx, var)
            %%% make indices and values for derivatives
            % find components with derivatives in var_name direction
            comp_in = find(obj.model.n_order > 0);
            n_comp = numel(comp_in);

            % number of parameters, when not taylor expanded
            npars_ne = obj.model.npars_not_expanded;
            assert(all(npars_ne == npars_ne(1)),...
                'Number of non-expanded parameters per component should match')
            npars_ne = npars_ne(1);

            % total number of cells
            n_cells = obj.mesh.ncells;

            % component names
            comp_names = obj.model.component_names;

            % holds indices of parameters to be used later. For each
            % derivative to var_name we will use the derivative of the 
            % parameter and parameter itself (order zero)
            f_par = nan(n_comp, npars_ne*n_cells, 2);
            for cc=1:numel(comp_in) % for each component
                cur_comp = comp_names{comp_in(cc)}; % name of component
                f_par(cc,:,1) = find(obj.find_par(1,cur_comp, var_name)); % dp/dvar
                f_par(cc,:,2) = find(obj.find_par(0,cur_comp)); % p0
            end

            % reshape to be able to easily search for neigbors in next step
            f_par = reshape(f_par,...
                n_comp,... % number of components
                npars_ne,... % number of not expanded parameters
                n_cells, ... % number of cells
                2); % derivatives of paramters to n and parameter itself
            
            % initialize column indices that will refer to derivative of
            % parameter and the parameter in the cell to the right and to
            % the left
            nb = obj.neighbors;
            has_neighbors = all(isfinite(nb(nb_idx,:)),1);
            n_in = sum(has_neighbors);
            col = nan(n_comp, npars_ne, n_in, 3);
            for cc = 1:numel(comp_in)
                col(cc, :, :, 1) = ...
                    f_par(cc, :, has_neighbors, 1); %dp/dn
                col(cc, :, :, 2) = ...
                    f_par(cc, :, nb(nb_idx(1), has_neighbors), 2); %p0: right
                col(cc, :, :, 3) = ...
                    f_par(cc, :, nb(nb_idx(2), has_neighbors), 2); %p0: left
            end

            % make row indices. Each derivative for each component is a
            % different equation with three terms
            
            row = repmat((1 : n_comp * npars_ne * n_in)', ...
            [1, 3]);

            % reshape row indices
            col = reshape(col,[],3);

            % make coefficients. Note that the coefficients for the
            % derivatives have a minus to the right and a plus to the lift,
            % because we want to negate the finite difference. Coefficients
            % are vary only with cells, and not with component or parameter
            dvar = var(nb(nb_idx(1), has_neighbors))-...
                 var(nb(nb_idx(2), has_neighbors));
            val = cat(3,...
                ones(1,n_in),...
                -1./dvar,...
                 1./dvar);

            % replicate and reshape values
            val = shiftdim(val, -1);
            val = repmat(val,[n_comp, npars_ne]);
            val = reshape(val, [], 3);
        end
    end
end