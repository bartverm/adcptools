classdef Regularization <...
        handle &... % handle object
        matlab.mixin.Heterogeneous &... % allow arrays of different subclasses
        helpers.ArraySupport % add array functionality

    properties
        par_names_tot (:,1) cell
        neighbors (:,1) cell
        domains (:,1) cell
        zb_0 (:,1) double
        zb_x (:,1) double
        zb_s (:,1) double
        zb_n (:,1) double
        C (1,:) cell = {};
    end

    methods(Sealed)
        function C = assemble_matrices(obj, solver)
            if ~isscalar(obj)
                C = obj.run_method('assemble_matrices', solver);
            else
                C = obj.assemble_matrix_private();
            end
        end
    end

    methods(Abstract, Access = protected)
        assemble_matrix_private(obj)
    end

    methods
        function obj = Regularization(bathy, xs, mesh, data_model, water_level)
            obj.par_names_tot = cell([mesh.ncells*sum(data_model.npars),1]);
            obj.neighbors = zeros([4, mesh.ncells]);
            obj.domains = zeros([1, mesh.ncells]);

            dx = 2; 
            dy = 2; % meters: for estimating bathymetric gradients using centered differences
            
            ds = 2;
            dn = 2;

            % Here, the specific orientation of cross-sections becomes
            % important.

            svec = xs.direction_orthogonal;
            nvec = xs.direction;

            for idx = 1:mesh.ncells % loop trough every cell
                for par = 1:sum(data_model.npars) % loop trough parameters within cell
                    obj.par_names_tot{idx, 1} = sprintf('cell %i: %s', idx, data_model.names{par});
                end
                [obj.neighbors(:,idx), obj.domains(1,idx)] = mesh.get_neighbors(idx);
                % Bathymetry and its gradients in along-channel and
                % cross-channel directions.
                obj.zb_0(1,idx) = mesh.zb_middle(mesh.col_to_cell(idx));

                % Centered difference in Cartesian x,y coordinates
                obj.zb_x(1,idx) = 1/(2*dx)*(bathy.get_depth(mesh.x_middle(mesh.col_to_cell(idx))+dx,mesh.y_middle(mesh.col_to_cell(idx)),t)-...
                    bathy.get_depth(mesh.x_middle(mesh.col_to_cell(idx))-dx,mesh.y_middle(mesh.col_to_cell(idx)),mesh.time));
                obj.zb_y(1,idx) = 1/(2*dy)*(bathy.get_depth(mesh.x_middle(mesh.col_to_cell(idx)),mesh.y_middle(mesh.col_to_cell(idx))+dy,t)-...
                    bathy.get_depth(mesh.x_middle(mesh.col_to_cell(idx)),mesh.y_middle(mesh.col_to_cell(idx))-dy,mesh.time));

                % Directional difference in cross-sectional s,n coordinates
                obj.zb_s(1,idx) = 1/(2*ds)*(bathy.get_depth(mesh.x_middle(mesh.col_to_cell(idx))+ds*svec(1), mesh.y_middle(mesh.col_to_cell(idx) + ds*svec(2)),t)-...
                    bathy.get_depth(mesh.x_middle(mesh.col_to_cell(idx)) - ds*svec(1), mesh.y_middle(mesh.col_to_cell(idx) - ds*svec(2)),mesh.time));
                obj.zb_n(1,idx) = 1/(2*dn)*(bathy.get_depth(mesh.x_middle(mesh.col_to_cell(idx))+dn*nvec(1),mesh.y_middle(mesh.col_to_cell(idx))+dn*nvec(2),t)-...
                    bathy.get_depth(mesh.x_middle(mesh.col_to_cell(idx)) - dn*nvec(1),mesh.y_middle(mesh.col_to_cell(idx)) - dn*svec(2),mesh.time));
            end
        end
    end
end