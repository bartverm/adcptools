classdef Regularization <...
        helpers.ArraySupport % add array functionality

% Regularization of solution of ADCP data on a mesh
%
%   This is a generic class to define regularizations to be used when
%   solving models of ADCP data on a mesh. 
%
%   regularization.Regularization properties:
%   bathy - Bathymetry object defining the bathymetry
%   xs - XSection object defining the cross-section
%   mesh - Mesh object defining the mesh on which data is to be solved
%   model - DataModel object defining the model to be solved with the data
%   weight - weight to be used for the regularization
%
%   regularization.Regularization read only properties:
%   C - regularization matrix
%   Cg - gramian of the regularization matrix, i.e. C'*C
%   rhs - right hand side term of the regularization
%   assembled - whether the above matrices are already assembled
%
%   regularization.Regularization method:
%   get_all_regs - returns all available regularizations
%   assemble_matrices - assemble the regularization matrices
%
%   see also: Solver, DataModel, regularization.Velocity,
%   regularization.TaylorBased

    properties(SetObservable)
        % Regularization/bathy
        %
        %   Bathymetry object defining the bathymetry around the mesh on
        %   which model is solved. Used to compute gradients in bathymetry
        %
        %   see also: Regularization
        bathy (1,1) Bathymetry = BathymetryScatteredPoints

        % Regularization/xs
        %
        %   XSection object defining the cross-section where mesh is
        %   defined on which model is being solved
        %
        %   see also: Regularization
        xs (1,1) XSection = XSection

        % Regularization/mesh
        %
        %   Mesh on which model is being solved
        %
        %   see also: Regularization
        mesh (1,1) Mesh = SigmaZetaMesh

        % Regularization/model
        %
        %   Data model that is to be solved
        %
        %   see also: Regularization
        model (1,1) DataModel = VelocityModel

        % Regularization/weight
        %
        %   Weights to use in the regularization. This can be a vector
        %   causing the Solver to produce different solutions. Setting
        %   weight to zero will result in disabling the regularization.
        %   The larger this value, the stronger the regularization is
        %   imposed on the solution.
        %
        %   see also: Regularization
        weight (1,:) double {mustBeNonempty,...
            mustBeFinite,...
            mustBeNonnegative} = 1
    end
    properties(SetAccess = protected)
        % Regularization/C
        %
        %   regularization matrix. Matrix is assembled on the fly when
        %   needed.
        %
        %   see also: Regularization, assembled, Cg, rhs
        C (:,:) double {helpers.mustBeSparse} = sparse(0)

        % Regularization/Cg
        %
        %  Gramian of the regularization matrix, i.e. C'*C.
        %  Matrix is assembled on the fly when needed.
        %
        %   see also: Regularization, assembled, C, rhs
        Cg (:,:) double {helpers.mustBeSparse} = sparse(0)

        % Regularization/rhs
        %
        %  Right hand side values for the regularization
        %  Matrix is assembled on the fly when needed.
        %
        %   see also: Regularization, assembled, Cg, C
        rhs (:,1) double {helpers.mustBeSparse} = sparse(0)

        % Regularization/assembled
        %
        %  logical value indicating whether matrices where already
        %  assembled.
        %  Matrices are only assembled when this value is false.
        %
        %   see also: Regularization, C, Cg, rhs
        assembled (1,1) logical = false
    end

    properties(GetAccess = protected, SetAccess = private, Dependent)
        Cg_weight (:,:) double
        rhs_weight (:,:) double
        neighbors (:,1) cell
        domains (:,1) cell
        zb0 (1,:) double
        zbxy (2,:) double
        zbsn (2,:) double
    end
    methods
        function val = get.C(obj)
            obj.assemble_matrices
            val = obj.C;
        end
        function val = get.Cg(obj)
            obj.assemble_matrices
            val = obj.Cg;
        end
        function val = get.rhs(obj)
            obj.assemble_matrices
            val = obj.rhs;
        end
        function neighbors = get.neighbors(obj)
            neighbors = obj.mesh.neighbors;
        end
        function domains = get.domains(obj)
            domains = obj.mesh.domains;
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
        function check_regpars(obj)
            reg_pars = {obj.weight};
            siz_reg = cellfun(@size,reg_pars,'UniformOutput',false);
            assert(isscalar(reg_pars) || isequal(siz_reg{:}),...
                'Weights of regression objects must have the same size')
        end
        function assemble_matrix_private(obj)
            obj.C = sparse(0);
        end
    end
    methods(Static)
        function obj = get_all_regs(varargin)
        % Return all available regularizations.
        %
        %   this returns Regularization, which is the default
        %   regularization. Object are different for e.g.
        %   regularizaton.Velocity, that will return all available velocity
        %   regularizations.
        %
        %   See also: Regularization, regularization.Velocity.get_all_regs
            obj = regularization.Regularization(varargin{:});
        end
    end
    methods(Static, Access = protected)
        function regs = reorganize_regs(regs)
            siz = cellfun(@size, regs, 'UniformOutput',false);
            assert(isequal(siz{:}),...
                'Size of generated regularizations should be equal')
            siz = siz{1};
            if isequal(siz, [1 1])
                regs = [regs{:}];
                return
            end
            out = cell(siz);
            for co = 1:numel(out)
                tmp = cellfun(@(x) x(co), regs, 'UniformOutput',false);
                tmp = [tmp{:}];
                out{co} = tmp;
            end
            regs = out;

        end
    end

    methods(Sealed)
        function assemble_matrices(obj)
        % assemble regularization matrices
        %
        %   obj.assemble_matrices(obj) constructs the regularization
        %   matrices if they have not been constructed yet, i.e. if
        %   obj.assembled is false. After successfull assembly, the
        %   assembled property is set to true;
        %
        %   see also: Regularization, assembled
            if ~isscalar(obj)
                obj.run_method('assemble_matrices');
            else
                if obj.assembled == false
                    % set assembled to true before assembling to avoid
                    % infinite recursion.
                    obj.assembled = true;
                    % try, catch, to makes sure assembled is set to false
                    % if assembly fails
                    name = class(obj);
                    try
                        obj.assemble_matrix_private(); 
                        obj.gramian_matrix();
                        disp(strcat("Assembled regularization matrix ", name(16:end)))
                    catch err
                        %obj.assembled = false;
                        %rethrow(err)
                        warning(err.message)
                        disp(strcat("Failed to assemble regularization matrix ", name(16:end), "-> using zero matrix"))
                        obj.C = sparse(0);
                        obj.Cg = sparse(0);
                        obj.rhs = sparse(0);
                    end
                end
            end
        end
    end
    methods(Access = protected)
        function gramian_matrix(obj)
            obj.Cg = obj.C'*obj.C;
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
    end
end
