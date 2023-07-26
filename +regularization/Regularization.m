classdef Regularization <...
        helpers.ArraySupport% &... % add array functionality
        % matlab.mixin.Heterogeneous  % allow arrays of different subclasses

    properties(SetObservable)
        bathy (1,1) Bathymetry = BathymetryScatteredPoints
        xs (1,1) XSection = XSection
        mesh (1,1) Mesh = SigmaZetaMesh
        model (1,1) DataModel = VelocityModel
        weight (1,:) double {mustBeNonempty,...
            mustBeFinite,...
            mustBeNonnegative} = 1

    end
    properties(SetAccess = protected)
        % Regularization/matrix
        C (:,:) double {helpers.mustBeSparse} = sparse(0)
        Cg (:,:) double {helpers.mustBeSparse} = sparse(0)
        rhs (:,1) double {helpers.mustBeSparse} = sparse(0)
        assembled (1,1) logical = false
    end

    properties(GetAccess = public, SetAccess = private, Dependent)
        Cg_weight (:,:) double
        rhs_weight (:,:) double
        names_all (1,:) cell
        neighbors (:,1) cell
        domains (:,1) cell
        zb0 (1,:) double
        zbxy (2,:) double
        zbsn (2,:) double
    end
    methods
        function val = get.Cg_weight(obj)
            if isscalar(obj)
                val = obj.Cg .* shiftdim(obj.weight, -1);
                return
            end
            obj.check_regpars;
            val = obj(1).Cg_weight;
            for co = 2:numel(obj)
                val = val + obj(co).Cg_weight;
            end
        end
    end
    methods(Access = protected)
        function check_regpars(obj)
            reg_pars = {obj.weight};
            siz_reg = cellfun(@size,reg_pars,'UniformOutput',false);
            assert(isscalar(reg_pars) || isequal(siz_reg{:}),...
                'Weights of regression objects must have the same size')
        end
        function assemble_matrix_private(obj)
            obj.assembled = true;
        end
    end
    methods(Static, Abstract)
        % Return all availabl regularizations
        get_all_regs()
    end

    methods(Sealed)
        function assemble_matrices(obj)
            if ~isscalar(obj)
                obj.run_method('assemble_matrices');
            else
                if obj.assembled == false
                    obj.assemble_matrix_private();
                    obj.gramian_matrix;
                    obj.assembled = true;
                end
            end
        end
    end
    methods
        function names_all = get.names_all(obj)
            % -> TaylorBasedRegular.
            flat_names = obj.flatten_names();
            cells_vec = cell([1,obj.mesh.ncells]);
            for idx = 1:obj.mesh.ncells % loop trough every cell
                cells_vec{1,idx} = sprintf('cell %i: ', idx);
            end
            names_all = helpers.kron_modified_cell(cells_vec, flat_names);
        end

        function neighbors = get.neighbors(obj)
            % get neighbors for each mesh cell -> Mesh
            neighbors = zeros([4,obj.mesh.ncells]);
            for idx = 1:obj.mesh.ncells % loop trough every cell
                [neighbors(:,idx), ~] = obj.mesh.get_neighbors(idx);
            end
        end

        function domains = get.domains(obj)
            % get domain -> Mesh
            domains = zeros([1,obj.mesh.ncells]);
            for idx = 1:obj.mesh.ncells % loop trough every cell
                [~, domains(1,idx)] = obj.mesh.get_neighbors(idx);
            end
        end

        function zb0 = get.zb0(obj)
            % get bed elevation
            zb0(1,:) = obj.mesh.zb_middle(obj.mesh.col_to_cell);
        end
        function zbxy = get.zbxy(obj)
            % compute dzb/dx dzb/dy for each mesh cell -> Mesh
            dx = 2;
            dy = 2; % meters: for estimating bathymetric gradients using centered differences
            % Centered difference in Cartesian x,y coordinates
            zbxy(1,:) = 1/(2*dx)*(obj.bathy.get_depth(obj.mesh.x_middle(obj.mesh.col_to_cell)+dx,obj.mesh.y_middle(obj.mesh.col_to_cell),obj.mesh.time)-...
                obj.bathy.get_depth(obj.mesh.x_middle(obj.mesh.col_to_cell)-dx,obj.mesh.y_middle(obj.mesh.col_to_cell),obj.mesh.time));
            zbxy(2,:) = 1/(2*dy)*(obj.bathy.get_depth(obj.mesh.x_middle(obj.mesh.col_to_cell),obj.mesh.y_middle(obj.mesh.col_to_cell)+dy,obj.mesh.time)-...
                obj.bathy.get_depth(obj.mesh.x_middle(obj.mesh.col_to_cell),obj.mesh.y_middle(obj.mesh.col_to_cell)-dy,obj.mesh.time));
        end
        function zbsn = get.zbsn(obj)
            % same as above, rotated to xs directions
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
        function IM = assemble_incidence(obj)
            % unused -> remove?
            obj.mesh = obj.mesh;
            IM = zeros(obj.mesh.ncells);
            for j = 1:obj.mesh.ncells % loop trough every cell
                nbreal = obj.neighbors(j,~isnan(obj.neighbors(j,:)));
                nbreal = nbreal(nbreal>j);
                IM(j, nbreal) = 1;
            end
            IM = IM + IM';
        end


        function gramian_matrix(obj)
            obj.Cg = obj.C'*obj.C;
        end
        
        function keep_idx = dom2keep_idx(obj)
            % -> ConsistencyRegularization
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
            flat_names = obj.model.all_names;
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
