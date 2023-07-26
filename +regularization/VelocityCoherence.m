classdef VelocityCoherence < regularization.Coherence &...
        regularization.Velocity
    methods(Access = protected)
        function assemble_matrix_private(obj)
            obj.assemble_matrix_private@regularization.Velocity
            obj.assemble_matrix_private@regularization.Coherence
        end
        function W = assemble_weights(obj)
            W = obj.assemble_weights@regularization.Coherence;

            % Extract diagonal
            w = spdiags(W,0);

            vertscale = max(obj.mesh.z_patch,[], 'all') - min(obj.mesh.z_patch, [], 'all');
            horscale = max(obj.mesh.n_patch, [], 'all') - min(obj.mesh.n_patch, [], 'all'); %Typical scales
            
            f_vert_vel = obj.find_par([],'w');

            w(f_vert_vel) = w(f_vert_vel)*horscale/vertscale;

            % Replace diagonal with new values
            W = spdiags(w, 0, W);
        end
    end
end