classdef Coherence < regularization.Regularization
    methods(Access = protected)
        function assemble_matrix_private(obj)
            Np = sum(obj.model.npars);
            ncells = obj.mesh.ncells;
            
            % index of four neighbors for each cell. If no neighbor is
            % present nb holds a NaN
            nb = obj.neighbors; 

            % number of neighbors for each cell
            num_neighbors = sum(isfinite(nb),1);

            % if there is no neighbor, make a face cell (at ncells+1) be
            % the neighbor, we'll remove these later, but needed to
            % vectorize computation
            nb(isnan(nb)) = ncells+1;
            
            % matrix holding the parameter numbers in the final parameter
            % vector for each cell. The matrix is Np x ncells + 1, which
            % includes the fake cell
            fpars = reshape(1:Np*(ncells+1), Np, ncells+1);

            % we make a vector that contains for each parameter, the index
            % of the four neighboring paramaters.
            col = cat(3,...
                fpars(:,nb(1,:)),...
                fpars(:,nb(2,:)),...
                fpars(:,nb(3,:)),...
                fpars(:,nb(4,:)));

            % vector containing the index of the parameter for which we are
            % constructing the equation, so in 3rd dimension it is repeated
            row = repmat(fpars(:,1:end-1), [1 1 4]);

            % value is 1/number_neighbors. This is the real number of
            % neighbors, excluding the face cell
            val = -ones(size(row))./num_neighbors;

            % we now remove the coefficients when the index refers to the
            % fake cell
            fbad = col > Np*ncells; % bad coefficients
            col = col(~fbad); 
            row = row(~fbad);
            val = val(~fbad);

            % construct the matrix
            obj.C = speye(Np*ncells) +... % identity matrix
                sparse(row, col, val, ncells*Np, ncells*Np);

            W = obj.assemble_weights();

            obj.C = W*obj.C;

        end
    end
end



