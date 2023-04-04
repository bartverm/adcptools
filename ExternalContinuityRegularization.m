classdef ExternalContinuityRegularization < Regularization
methods (Access = protected)
        function C2 = assemble_matrix_private(obj)

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

end
end