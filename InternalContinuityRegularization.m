classdef InternalContinuityRegularization < Regularization
    methods(Access=protected)
        function C1=assemble_matrix_private(obj)
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
    end
end