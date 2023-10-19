classdef VelocityCoherence < regularization.Coherence &...
        regularization.Velocity
% Weakly imposes spatial coherence between parameters in a velocity model
%
% see also: Coherence, Regularization
    methods(Access = protected)
        function W = assemble_weights(obj)
            W = obj.assemble_weights@regularization.Coherence;

            % Extract diagonal
            w = spdiags(W,0);

            vertscale = max(obj.mesh.z_patch,[], 'all') - min(obj.mesh.z_patch, [], 'all');
            horscale = max(obj.mesh.n_patch, [], 'all') - min(obj.mesh.n_patch, [], 'all'); %Typical scales
            if isa(obj.model, 'TaylorModel') || isa(obj.model, 'TidalModel')
                f_vert_vel = obj.find_par(component = 'w');
            else
                f_vert_vel = 3:3:size(w,1);
            end
            w(f_vert_vel) = w(f_vert_vel)*horscale/vertscale;

            % Replace diagonal with new values
            W = spdiags(w, 0, W);
        end
    end
end