classdef Coherence < regularization.Regularization &...
        regularization.TaylorBased
% Impose spatial coherence between neighboring cells
%
%   This regularization imposes spatial coherence of model parameters
%   between the cells. A larger weight will impose a larger spatial
%   coherence resulting in a smoother solution in space.
%
%   see also: VelocityCoherence, Regularization
    methods(Access = protected)
        function assemble_matrix_private(obj)
            Np = sum(obj.model.npars);
            ncells = obj.mesh.ncells;
            
            % index of four neighbors for each cell. If no neighbor is
            % present nb holds a NaN
            nb = obj.neighbors'; 

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
            C = speye(Np*ncells) +... % identity matrix
                sparse(row, col, val, ncells*Np, ncells*Np);

            W = obj.assemble_weights();

            obj.C = W*C;

            obj.rhs = sparse([],[],[],ncells*Np,1,0);

        end

        function W = assemble_weights(obj)
            % scaling factors based on characteristic scales

            npars = sum(obj.model.npars)*obj.mesh.ncells;
            w = ones([npars,1]);

            if ~isa(obj.model,'TaylorModel')
                W = spdiags(w, 0, npars, npars);
                return
            end

            % Apply enhanced regularization for small variables using
            % characteristic spatial scales

            horscale = max(obj.mesh.n_patch, [], 'all') - min(obj.mesh.n_patch, [], 'all'); %Typical scales

            % Automatically assign weights to smoothness of different parameters
            % Important due to orientation of main flow
            f_horz_der = obj.find_par(order = 1, variable = {'n','s'});

            w(f_horz_der) = w(f_horz_der)*horscale;
            W = spdiags(w, 0, npars, npars);
        end
    end
    methods(Access = protected)
        function val = get_min_order(obj)
            nc = obj.model.ncomponents;
            val = zeros(5,nc);
        end
    end
end



