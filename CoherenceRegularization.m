classdef CoherenceRegularization < Regularization
    methods(Access = protected)
        function C3 = assemble_matrix_private(obj)

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

            W = obj.assemble_weights();

            C3 = W*C3;

        end
    end
end



