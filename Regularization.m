classdef Regularization <...
        handle &... % handle object
        helpers.ArraySupport % add array functionality
    %matlab.mixin.Heterogeneous &... % allow arrays of different subclasses

    properties
        bathy
        xs (1,1) XSection
        mesh (1,1) SigmaZetaMesh
        model (1,1) DataModel
        C (1,:) cell = {sparse(0), sparse(0), sparse(0), sparse(0), sparse(0)}
        Cg (1,:) cell = {sparse(0), sparse(0), sparse(0), sparse(0), sparse(0)}
        rhs (:,1) double = sparse(0);
        assembled (1,1) logical = false
    end

    properties(Dependent)
        names_all (1,:) cell
        neighbors (:,1) cell
        domains (:,1) cell

        zb0 (1,:) double
        zbxy (2,:) double
        zbsn (2,:) double

    end

    methods
        function obj = Regularization(varargin)
            % Constructor
            for ia = 1:2:nargin
                obj.(varargin{ia}) = varargin{ia+1};
            end
        end

        function assemble_matrices(obj)
            if(isa(obj.model,"TaylorModel"))
                if (obj.model.s_order(1) > 0) && (obj.model.n_order(2) > 0)...
                        && (obj.model.sigma_order(3) > 0 || obj.model.z_order(3) > 0)
                    obj.C{1} = obj.assemble_continuity_internal();
                else
                    warning('No internal continuity matrix assembled: Include higher order Taylor expansion')
                    obj.C{1} = sparse(0);
                end
                if (obj.model.s_order(1) > 0)
                    obj.C{2} = obj.assemble_continuity_external();
                else
                    warning('No external continuity matrix assembled: Include alongchannel Taylor expansion')
                    obj.C{2} = sparse(0);
                end
                if (all(obj.model.n_order > 0)) && (all(obj.model.sigma_order > 0) || all(obj.model.z_order > 0))
                    % This condition may be relaxed a bit.
                    obj.C{4} = obj.assemble_consistency();
                else
                    warning('No consistency matrix assembled: Fit a sufficient number of Taylor terms')
                    obj.C{4} = sparse(0);
                end
                if (all(obj.model.sigma_order > 0) || all(obj.model.z_order > 0))
                    [obj.C{5}, obj.rhs] = obj.assemble_kinematic();
                else
                    warning('No kinematic boundary condition matrix assembled: Include higher order Taylor expansion in sigma/z direction')
                    obj.C{5} = sparse(0);
                    obj.rhs=sparse(0);
                end
            else
                warning("To impose constraints on continuity and boundary conditions, a Taylor model is required. Only assembling coherence operator.");
                obj.C{1} = sparse(0);
                obj.C{2} = sparse(0);
                obj.C{4} = sparse(0);
                obj.C{5} = sparse(0);
                obj.rhs = sparse(0);
            end

            obj.C{3} = obj.assemble_coherence();
            obj.gramian_matrices()
            obj.assembled = true;
        end

        function names_all = get.names_all(obj)
            flat_names = obj.flatten_names();
            cells_vec = cell([1,obj.mesh.ncells]);
            for idx = 1:obj.mesh.ncells % loop trough every cell
                cells_vec{1,idx} = sprintf('cell %i: ', idx);
            end
            names_all = helpers.kron_modified_cell(cells_vec, flat_names);
        end

        function neighbors = get.neighbors(obj)
            neighbors = zeros([4,obj.mesh.ncells]);
            for idx = 1:obj.mesh.ncells % loop trough every cell
                [neighbors(:,idx), ~] = obj.mesh.get_neighbors(idx);
            end
        end

        function domains = get.domains(obj)
            domains = zeros([1,obj.mesh.ncells]);
            for idx = 1:obj.mesh.ncells % loop trough every cell
                [~, domains(1,idx)] = obj.mesh.get_neighbors(idx);
            end
        end

        function zb0 = get.zb0(obj)
            zb0(1,:) = obj.mesh.zb_middle(obj.mesh.col_to_cell);
        end
        function zbxy = get.zbxy(obj)
            dx = 2;
            dy = 2; % meters: for estimating bathymetric gradients using centered differences
            % Centered difference in Cartesian x,y coordinates
            zbxy(1,:) = 1/(2*dx)*(obj.bathy.get_depth(obj.mesh.x_middle(obj.mesh.col_to_cell)+dx,obj.mesh.y_middle(obj.mesh.col_to_cell),obj.mesh.time)-...
                obj.bathy.get_depth(obj.mesh.x_middle(obj.mesh.col_to_cell)-dx,obj.mesh.y_middle(obj.mesh.col_to_cell),obj.mesh.time));
            zbxy(2,:) = 1/(2*dy)*(obj.bathy.get_depth(obj.mesh.x_middle(obj.mesh.col_to_cell),obj.mesh.y_middle(obj.mesh.col_to_cell)+dy,obj.mesh.time)-...
                obj.bathy.get_depth(obj.mesh.x_middle(obj.mesh.col_to_cell),obj.mesh.y_middle(obj.mesh.col_to_cell)-dy,obj.mesh.time));
        end
        function zbsn = get.zbsn(obj)

            ds = 2;
            dn = 2;

            % Here, the specific orientation of cross-sections becomes
            % important.

            svec = obj.xs.direction_orthogonal;
            nvec = obj.xs.direction;
            % Directional difference in cross-sectional s,n coordinates
            zbsn(1,:) = 1/(2*ds)*(obj.bathy.get_depth(obj.mesh.x_middle(obj.mesh.col_to_cell)+ds*svec(1), obj.mesh.y_middle(obj.mesh.col_to_cell) + ds*svec(2),obj.mesh.time)-...
                obj.bathy.get_depth(obj.mesh.x_middle(obj.mesh.col_to_cell) - ds*svec(1), obj.mesh.y_middle(obj.mesh.col_to_cell) - ds*svec(2),obj.mesh.time));
            zbsn(2,:) = 1/(2*dn)*(obj.bathy.get_depth(obj.mesh.x_middle(obj.mesh.col_to_cell)+dn*nvec(1),obj.mesh.y_middle(obj.mesh.col_to_cell)+dn*nvec(2),obj.mesh.time)-...
                obj.bathy.get_depth(obj.mesh.x_middle(obj.mesh.col_to_cell) - dn*nvec(1),obj.mesh.y_middle(obj.mesh.col_to_cell) - dn*svec(2),obj.mesh.time));
        end


    end
    methods(Access = protected)
        function C1 = assemble_continuity_internal(obj)
            % Function that assembles cell-based continuity equation
            wl = obj.bathy.water_level;

            nscale = 1; % B and L are possible lateral and longitudinal scaling factors
            sscale = 1;

            D0 = obj.get_subtidal_depth();
            const_names = obj.get_const_names(); % Cell array
            D0s = -obj.zbsn(1,:);
            D0n = -obj.zbsn(2,:);
            sig = obj.mesh.sig_center;
            Cj = cell([obj.mesh.ncells,1]);
            % C1 is block diagonal and can thus be assembled per cell
            pnames = obj.flatten_names();

            col = cell([1, numel(const_names)]);
            term = cell([1, numel(const_names)]);
            for eq = 1:numel(const_names)
                col{eq} = [find(strcmp(pnames, ['d^1u/dx^1', const_names{eq}])) , ...
                    find(strcmp(pnames, ['d^1u/dsig^1', const_names{eq}])) , ...
                    find(strcmp(pnames, ['d^1v/dy^1', const_names{eq}])) , ...
                    find(strcmp(pnames, ['d^1v/dsig^1', const_names{eq}])) , ...
                    find(strcmp(pnames, ['d^1w/dsig^1', const_names{eq}]))];
                if eq > 1 % For tidal constituents, correlations between water level and velocity in sigma coordinates have to be included.
                    col{eq}((numel(col{eq})+1):(numel(col{eq})+2)) = [find(strcmp(pnames, ['d^1u/dx^1', const_names{1}])),...
                        find(strcmp(pnames, ['d^1v/dy^1', const_names{1}]))];
                end
            end

            for cell_idx = 1:obj.mesh.ncells
                Cj{cell_idx} = zeros([numel(const_names), sum(obj.model.npars)]);
                for eq = 1:numel(const_names)
                    term{eq} = [nscale*D0(cell_idx), nscale*(1-sig(cell_idx))*D0s(cell_idx), sscale*D0(cell_idx), sscale*(1-sig(cell_idx))*D0n(cell_idx), sscale*nscale];
                    if eq > 1                    
                        term{eq}((numel(term{eq})+1):(numel(term{eq})+2)) = [nscale*wl.parameters(eq),...
                        nscale*wl.parameters(eq)];
                    end
                    Cj{cell_idx}(eq,col{eq}) = term{eq};
                end
            end
            C1 = helpers.spblkdiag(Cj{:});
        end

        function C2 = assemble_continuity_external(obj)

            wl = obj.bathy.water_level;
            D0 = obj.get_subtidal_depth();
            const_names = obj.get_const_names(); % Cell array
            D0s = -obj.zbsn(1,:);
            D0n = -obj.zbsn(2,:);
            
            sig_center = obj.mesh.sig_center;
            n_center = obj.mesh.n_middle(obj.mesh.col_to_cell);
            nscale = 1;
            sscale = 1;

            nb = obj.neighbors;
            dom = obj.domains;
            rows = []; cols = []; terms = [];
            pnames = obj.names_all;
            row_idx = 0;
            for cell_idx = 1:obj.mesh.ncells
                %cell_idx
                if dom(cell_idx) == 0
                    % Extended continuity only works on interior of domain.
                    dsig = sig_center(nb(2, cell_idx))-sig_center(nb(4, cell_idx));
                    dn = n_center(nb(1, cell_idx))-n_center(nb(3, cell_idx));
                    col = cell([1,numel(const_names)]);
                    term = cell([1,numel(const_names)]);
                    row = cell([1,numel(const_names)]);
                    for eq = 1:numel(const_names)
                        col{eq} = [find(strcmp(pnames, ['cell ',num2str(cell_idx),': d^1u/dx^1', const_names{eq}])) , ...
                            find(strcmp(pnames, ['cell ',num2str(nb(2, cell_idx)),': u0', const_names{eq}])) , ...
                            find(strcmp(pnames, ['cell ',num2str(nb(4, cell_idx)),': u0', const_names{eq}])) , ...
                            find(strcmp(pnames, ['cell ',num2str(nb(1, cell_idx)),': v0', const_names{eq}])) , ...
                            find(strcmp(pnames, ['cell ',num2str(nb(3, cell_idx)),': v0', const_names{eq}])) , ...
                            find(strcmp(pnames, ['cell ',num2str(nb(2, cell_idx)),': v0', const_names{eq}])) , ...
                            find(strcmp(pnames, ['cell ',num2str(nb(4, cell_idx)),': v0', const_names{eq}])) , ...
                            find(strcmp(pnames, ['cell ',num2str(nb(2, cell_idx)),': w0', const_names{eq}])) , ...
                            find(strcmp(pnames, ['cell ',num2str(nb(4, cell_idx)),': w0', const_names{eq}]))];

                        term{eq} = [nscale*D0(cell_idx), nscale*(1-sig_center(cell_idx))*D0s(cell_idx)/dsig, -nscale*(1-sig_center(cell_idx))*D0s(cell_idx)/dsig,...
                            sscale*D0(cell_idx)/dn, -sscale*D0(cell_idx)/dn, sscale*(1-sig_center(cell_idx))*D0n(cell_idx)/dsig, -sscale*(1-sig_center(cell_idx))*D0n(cell_idx)/dsig,...
                            sscale*nscale/dsig, -sscale*nscale/dsig];
                        if eq > 1 % For tidal constituents, correlations between water level and velocity in sigma coordinates have to be included.
                            col{eq}(numel(col{eq}):(numel(col{eq})+1)) = [find(strcmp(pnames, ['cell ',num2str(cell_idx),': d^1u/dx^1', const_names{1}])),...
                                find(strcmp(pnames, ['cell ',num2str(cell_idx),': d^1u/dx^1', const_names{1}]))];
                            term{eq}(numel(term{eq}):(numel(term{eq})+1)) = [nscale*wl.parameters(eq),...
                                nscale*wl.parameters(eq)];
                        end
                        row{eq} = row_idx + eq*ones(size(col{eq}));
                    end
                    rows = [rows [row{:}]];
                    cols = [cols [col{:}]];
                    terms = [terms [term{:}]];
                    row_idx = max(rows);
                end
            end
            C2 = sparse(rows, cols, terms, obj.mesh.ncells*numel(const_names), obj.mesh.ncells*sum(obj.model.npars));
        end

        function C3 = assemble_coherence(obj)

            Np = sum(obj.model.npars);
            Diag = speye(Np*obj.mesh.ncells);
            rows = []; cols = []; vals = [];
            nb = obj.neighbors;

            for cell_idx = 1:obj.mesh.ncells %rows
                %cell_idx
                nb_ = nb(:,cell_idx);
                nb_ = nb_(~isnan(nb_)); 
                nnb = length(nb_);
                for nb_idx = 1:nnb
                    col = ((nb_(nb_idx)-1)*Np+1):(nb_(nb_idx)*Np);
                    row = (cell_idx-1)*Np+1:cell_idx*Np;
                    val = -1/nnb*ones(1,Np);
                    rows = [rows row];
                    cols = [cols col];
                    vals = [vals val];
                end
            end
            C3 = Diag + sparse(rows, cols, vals, obj.mesh.ncells*Np, obj.mesh.ncells*Np);

            W = obj.assemble_weights;

            C3 = W*C3;

        end

        function C4 = assemble_consistency(obj)

            const_names = obj.get_const_names(); % Cell array

            nb = obj.neighbors;
            rows = []; cols = []; terms = [];
            pnames = obj.names_all;
            row_idx = 1;
            [dn, dsig] = obj.dom2dndsig();
            keep_idx = obj.dom2keep_idx();
            for cell_idx = 1:obj.mesh.ncells
                for eq = 1:numel(const_names)
                    col = [obj.findn(pnames, ['cell ',num2str(cell_idx),': d^1u/dy^1', const_names{eq}]) , ...
                        obj.findn(pnames, ['cell ',num2str(nb(1, cell_idx)),': u0', const_names{eq}]) , ...
                        obj.findn(pnames, ['cell ',num2str(nb(3, cell_idx)),': u0', const_names{eq}]) , ...
                        obj.findn(pnames, ['cell ',num2str(cell_idx),': d^1v/dy^1', const_names{eq}]) , ...
                        obj.findn(pnames, ['cell ',num2str(nb(1, cell_idx)),': v0', const_names{eq}]) , ...
                        obj.findn(pnames, ['cell ',num2str(nb(3, cell_idx)),': v0', const_names{eq}]) , ...
                        obj.findn(pnames, ['cell ',num2str(cell_idx),': d^1w/dy^1', const_names{eq}]) , ...
                        obj.findn(pnames, ['cell ',num2str(nb(1, cell_idx)),': w0', const_names{eq}]) , ...
                        obj.findn(pnames, ['cell ',num2str(nb(3, cell_idx)),': w0', const_names{eq}]) , ...
                        obj.findn(pnames, ['cell ',num2str(cell_idx),': d^1u/dsig^1', const_names{eq}]) , ...
                        obj.findn(pnames, ['cell ',num2str(nb(2, cell_idx)),': u0', const_names{eq}]) , ...
                        obj.findn(pnames, ['cell ',num2str(nb(4, cell_idx)),': u0', const_names{eq}]) , ...
                        obj.findn(pnames, ['cell ',num2str(cell_idx),': d^1v/dsig^1', const_names{eq}]) , ...
                        obj.findn(pnames, ['cell ',num2str(nb(2, cell_idx)),': v0', const_names{eq}]) , ...
                        obj.findn(pnames, ['cell ',num2str(nb(4, cell_idx)),': v0', const_names{eq}]) , ...
                        obj.findn(pnames, ['cell ',num2str(cell_idx),': d^1w/dsig^1', const_names{eq}]) , ...                        
                        obj.findn(pnames, ['cell ',num2str(nb(2, cell_idx)),': w0', const_names{eq}]), ...
                        obj.findn(pnames, ['cell ',num2str(nb(4, cell_idx)),': w0', const_names{eq}])];

                    row = row_idx + [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5]; % Row incremental index, always the same
                    term = [1, -1/(dn(cell_idx)), 1/(dn(cell_idx)),...
                        1, -1/(dn(cell_idx)), 1/(dn(cell_idx)),...
                        1, -1/(dn(cell_idx)), 1/(dn(cell_idx)),...
                        1, -1/(dsig(cell_idx)), 1/(dsig(cell_idx)),...
                        1, -1/(dsig(cell_idx)), 1/(dsig(cell_idx)),...
                        1, -1/(dsig(cell_idx)), 1/(dsig(cell_idx))]; % Value, always the same (but look at order of magnitudes difference between the values...)
                    
                    rows = [rows row(keep_idx{cell_idx})];
                    cols = [cols col(keep_idx{cell_idx})];%(keep_idx{cell_idx})];
                    terms = [terms term(keep_idx{cell_idx})];
                    if ~isempty(rows)
                        row_idx = max(rows) + 1;
                    end
                end
            end
            C4 = sparse(rows, cols, terms, 6*obj.mesh.ncells*numel(const_names), obj.mesh.ncells*sum(obj.model.npars));
        end

        function [C5, rhsvec] = assemble_kinematic(obj)

            % Function that assembles cell-based kinematic boundary conditions.

            const_names = obj.get_const_names(); % Cell array
            zbs = obj.zbsn(1,:);
            zbn = obj.zbsn(2,:);
            dom = obj.domains;
            col = cell([1, numel(const_names)]);
            pnames = obj.flatten_names();
            for eq = 1:numel(const_names)
                col{eq} = [find(strcmp(pnames, ['u0', const_names{eq}])) , ...
                    find(strcmp(pnames, ['d^1u/dsig^1', const_names{eq}])) , ...
                    find(strcmp(pnames, ['v0', const_names{eq}])) , ...
                    find(strcmp(pnames, ['d^1v/dsig^1', const_names{eq}])) , ...
                    find(strcmp(pnames, ['w0', const_names{eq}])) , ...
                    find(strcmp(pnames, ['d^1w/dsig^1', const_names{eq}]))];
            end
            % Now, the indices of the relevant terms are known. Proceed to fill in
            Cj = cell([obj.mesh.ncells,1]);
            rhsj = cell([obj.mesh.ncells,1]);

            for cell_idx = 1:obj.mesh.ncells
                Cj{cell_idx} = zeros([numel(const_names),sum(obj.model.npars)]);
                rhsj{cell_idx} = zeros([numel(const_names),1]);

                if dom(cell_idx)==0 || dom(cell_idx)==1 || dom(cell_idx) == 5 % Internal cells (in terms of sigma)
                    % Do nothing -> No lateral boundary conditions
                else
                    bot = (dom(cell_idx) ==6 || dom(cell_idx)==7|| dom(cell_idx)==8);  %Bottom cells
                    surf = (dom(cell_idx) ==2 || dom(cell_idx)==3|| dom(cell_idx)==4); %Surface cells
                    if surf
                        dsig = 1 - obj.mesh.sig_center(cell_idx); % per definition of sigma. Positive
                        terms = [0, 0, 0, 0, 1, dsig];
                        rhsj{cell_idx} = obj.water_level2rhs_element;
                    elseif bot % Bottom cells
                        dsig = obj.mesh.sig_center(cell_idx); %
                        terms = [zbs(cell_idx), -dsig*zbs(cell_idx), zbn(cell_idx), -dsig*zbn(cell_idx), -1, dsig];
                    end
                    for eq = 1:numel(const_names)
                        Cj{cell_idx}(eq, col{eq}) = terms;
                    end
                end
            end
            C5 = helpers.spblkdiag(Cj{:});
            rhsvec = cell2mat(rhsj);
        end

        function IM = assemble_incidence(obj)
            obj.mesh = obj.mesh;
            IM = zeros(obj.mesh.ncells);
            for j = 1:obj.mesh.ncells % loop trough every cell
                nbreal = obj.neighbors(j,~isnan(obj.neighbors(j,:)));
                nbreal = nbreal(nbreal>j);
                IM(j, nbreal) = 1;
            end
            IM = IM + IM';
        end

        function W = assemble_weights(obj)

            par_names = obj.flatten_names;
            Np = sum(obj.model.npars);
            w = ones([Np,1]);

            % Apply enhanced regularization for small singular value features using
            % characteristic spatial scales

            vertscale = max(obj.mesh.z_patch,[], 'all') - min(obj.mesh.z_patch, [], 'all');
            horscale = max(obj.mesh.n_patch, [], 'all') - min(obj.mesh.n_patch, [], 'all'); %Typical scales

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
            wj = cell(1,obj.mesh.ncells);
            wj(:)= {diag(w)};
            W = helpers.spblkdiag(wj{:});
        end

        function gramian_matrices(obj)
            obj.Cg = cell(size(obj.C));
            for idx = 1:length(obj.C)
                obj.Cg{idx} = obj.C{idx}'*obj.C{idx};
            end
        end
        
        function keep_idx = dom2keep_idx(obj)
            keep_idx = cell([obj.mesh.ncells,1]);
            dom = obj.domains;
            for cell_idx = 1:obj.mesh.ncells
                if dom(cell_idx) == 0
                    keep_idx{cell_idx} = 1:18;
                elseif dom(cell_idx) == 1
                    keep_idx{cell_idx} = 10:18;
                elseif dom(cell_idx) == 2
                    keep_idx{cell_idx} = [];
                elseif dom(cell_idx) == 3
                    keep_idx{cell_idx} = 1:9;
                elseif dom(cell_idx) == 4
                    keep_idx{cell_idx} = [];
                elseif dom(cell_idx) == 5
                    keep_idx{cell_idx} = 10:18;
                elseif dom(cell_idx) == 6
                    keep_idx{cell_idx} = [];
                elseif dom(cell_idx) == 7
                    keep_idx{cell_idx} = 1:9;
                elseif dom(cell_idx) == 8
                    keep_idx{cell_idx} = [];
                end
            end
        end

        function const_names = get_const_names(obj)
            if isa(obj.model, 'TidalModel')
                const_names{1} = ': M0'; % Subtidal always included
                for c = 1:numel(obj.model.constituents)
                    const_names{2*c} = [': ', obj.model.constituents{c}, 'a'];
                    const_names{2*c+1} = [': ', obj.model.constituents{c}, 'b'];
                end
            else
                const_names{1} = [];
            end

        end
        function D0 = get_subtidal_depth(obj)
            if isprop(obj.bathy.water_level, 'model') % If wl has a tidal model
                if isa(obj.bathy.water_level.model, 'TidalModel')
                    D0 = obj.bathy.water_level.parameters(1) - obj.zb0; % Subtidal depth
                end
            else
                D0 = obj.bathy.water_level.level - obj.zb0;
            end
        end

        function flat_names = flatten_names(obj)
            flat_names = [obj.model.names{:}];
        end

        function res = findn(obj, cell_of_str, str)
            res = find(strcmp(cell_of_str, str));
            if isempty(res)
                res = nan;
            end
        end

        function rhsj = water_level2rhs_element(obj)
            const_names = obj.get_const_names();
            wl = obj.bathy.water_level;
            
            rhsj = zeros([numel(const_names),1]);
            rhsj(1,1) = 0; % subtidal
            if numel(const_names) > 1
                omega = wl.model.get_omega;
                for eq = 1:numel(wl.model.constituents)
                    rhsj(2*eq, 1) = omega(eq)*wl.parameters(2*eq+1);
                    rhsj(2*eq + 1, 1) = -omega(eq)*wl.parameters(2*eq);
                end
            end
        end

        function [dn, dsig] = dom2dndsig(obj)
            nb = obj.neighbors;
            dom = obj.domains;
            dsig = zeros(size(dom));
            dn = zeros(size(dom));

            sig_center = obj.mesh.sig_center;
            n_center = obj.mesh.n_middle(obj.mesh.col_to_cell);

            for cell_idx = 1:obj.mesh.ncells
                if dom(cell_idx)==0 || dom(cell_idx)==1 || dom(cell_idx) == 5 % Internal cells (in terms of sigma)
                    dsig(cell_idx) = sig_center(nb(2, cell_idx))-sig_center(nb(4, cell_idx)); % central difference
                elseif dom(cell_idx) ==2 || dom(cell_idx) == 3|| dom(cell_idx)==4 %Surface cells (deprecated: replaced by boundary condition matrix)
                    dsig(cell_idx) = sig_center(cell_idx)-sig_center(nb(4, cell_idx)); % one-sided difference
                else                                                        % Bottom cells (deprecated: replaced by boundary condition matrix)
                    dsig(cell_idx) = sig_center(nb(2, cell_idx))-sig_center(cell_idx); % one-sided difference
                end

                if dom(cell_idx)==0 || dom(cell_idx)==3 || dom(cell_idx) == 7 % Internal cells (laterally)
                    dn(cell_idx) = n_center(nb(1, cell_idx))-n_center(nb(3, cell_idx)); % central difference
                elseif dom(cell_idx) ==1 || dom(cell_idx)==2|| dom(cell_idx)==8                           % Right side of domain (deprecated)
                    dn(cell_idx) = n_center(cell_idx)-n_center(nb(3, cell_idx)); % one-sided difference
                else                                                        % Left side of domain (deprecated)
                    dn(cell_idx) = n_center(nb(1, cell_idx))-n_center(cell_idx); % one-sided difference
                end
            end
        end
    end
end
