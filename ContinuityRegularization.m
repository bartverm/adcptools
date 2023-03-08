classdef ContinuityRegularization < Regularization
    properties
        C1 = 0;
        C2 = 0;
    end
    methods(Access = protected)
        function assemble_matrix_private(obj, solver)
            assert(isa(data_model,"TaylorModel"), "To constrain on continuity, a Taylor data_model is required");
            if (data_model.s_order(1) > 0) && (data_model.n_order(2) > 0)...
                    && (data_model.sigma_order(3) > 0 || data_model.z_order(3) > 0)
                obj.C1 = assemble_continuity_internal(obj, solver);
            else
                warning('No internal continuity matrix assembled: Include higher order Taylor expansion')
            end
            if (data_model.s_order(1) > 0)
                obj.C2 = assemble_continuity_external(obj, solver);
            else
                warning('No external continuity matrix assembled: Include alongchannel Taylor expansion')
            end
        end

        function C = assemble_continuity_internal(obj, bathy, xs, mesh, data_model, water_level)
            % Function that assembles cell-based continuity equation
            data_model = data_model;
            wl = water_level;

            nscale = 1; % B and L are possible lateral and longitudinal scaling factors
            sscale = 1;
            Cj = cell(mesh.ncells);
            for idx = 1:mesh.ncells
                sig = mesh.sig_center(idx);

                D0 = water_level.parameters(1) - mesh.zb_middle(mesh.col_to_cell(idx)); % Subtidal depth
                D0s = -obj.zb_s(idx);                                                                 % Subtidal depth gradient
                D0n = -obj.zb_n(idx);

                Cj{idx} = zeros([1+2*numel(wl.constituents), sum(data_model.npars)]);

                subtidal_idx = [find(strcmp(data_model.names, 'd^1u/dx^1: M0A')) , ...
                    find(strcmp(data_model.names, 'd^1u/dsig^1: M0A')) , ...
                    find(strcmp(data_model.names, 'd^1v/dy^1: M0A')) , ...
                    find(strcmp(data_model.names, 'd^1v/dsig^1: M0A')) , ...
                    find(strcmp(data_model.names, 'd^1w/dsig^1: M0A'))];

                subtidal_terms = [nscale*D0, nscale*(1-sig)*D0s, sscale*D0, sscale*(1-sig)*D0n, sscale*nscale];
                Cj{idx}(1,subtidal_idx) = subtidal_terms;

                % Second: Tidal equations for each constituent

                if isa(data_model,"TidalModel")
                    tidal_idx = cell([numel(data_model.constituentsU),2]);
                    tidal_terms = cell([numel(data_model.constituentsU),2]);
                    for tid_idx = 1:numel(data_model.constituentsU) %Fundamental choice: All constituents are the same!!
                        const = data_model.constituentsU{tid_idx};

                        tidal_idx{tid_idx,1} = [find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1u/dx^1: ', const,'A'))) , ...
                            find(strcmp(data_model.names, 'd^1u/dx^1: M0A')) , ...
                            find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1u/dsig^1: ', const,'A'))) ,...
                            find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1v/dy^1: ', const,'A'))) , ...
                            find(strcmp(data_model.names, 'd^1v/dy^1: M0A')) , ...
                            find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1v/dsig^1: ', const,'A'))) , ...
                            find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1w/dsig^1: ', const,'A')))];

                        tidal_terms{tid_idx,1} = [nscale*D0, nscale*wl.parameters(2*tid_idx), nscale*(1-sig)*D0s, sscale*D0, nscale*wl.parameters(2*tid_idx), sscale*(1-sig)*D0n, sscale*nscale];
                        Cj{idx}(2*tid_idx,tidal_idx{tid_idx,1}) = tidal_terms{tid_idx,1};

                        tidal_idx{tid_idx,2} = [find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1u/dx^1: ', const,'B'))) , ...
                            find(strcmp(data_model.names, 'd^1u/dx^1: M0A')) , ...
                            find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1u/dsig^1: ', const,'B'))) , ...
                            find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1v/dy^1: ', const,'B'))) , ...
                            find(strcmp(data_model.names, 'd^1v/dy^1: M0A')) , ...
                            find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1v/dsig^1: ', const,'B'))) , ...
                            find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1w/dsig^1: ', const,'B')))];

                        tidal_terms{tid_idx,2} = [nscale*D0, nscale*wl.parameters(2*tid_idx+1), nscale*(1-sig)*D0s, sscale*D0, nscale*wl.parameters(2*tid_idx+1), sscale*(1-sig)*D0n, sscale*nscale];
                        Cj{idx}(2*tid_idx+1,tidal_idx{tid_idx,2}) = tidal_terms{tid_idx,2};
                    end
                end
            end
            C = spblkdiag(Cj{:}); % Global matrices can be assembled from element matrices (C1 is block diagonal)
        end

        function C = assemble_continuity_external(obj, bathy, xs, mesh, data_model, water_level)

            mesh = mesh;
            data_model = data_model;
            wl = water_level;

            rows = []; cols = []; vals = [];
            row_idx = 1;
            for idx = 1:mesh.ncells
                for tid_idx = 1:length(wl.names)
                    tid_name = wl.names{tid_idx}(end-2:end); % Last three letters
                    [row, col, val] = obj.idx2element_mat(idx, tid_name, solver);
                    rows = [rows row_idx + row];
                    cols = [cols col];
                    vals = [vals val];
                    row_idx = max(rows) + 1;           
                end
            end
            C = sparse(rows, cols, vals, mesh.ncells*(2*length(data_model.constituentsU)+1), mesh.ncells*sum(data_model.npars));
        end

        function  [row, col, val] =  idx2element_mat(obj, idx, tid_name, solver)
            sig_center = mesh.sig_center;
            n_center = mesh.n_middle(mesh.col_to_cell);
            nscale = 1;
            sscale = 1;

            D0 = water_level.parameters(1) - mesh.zb_middle(mesh.col_to_cell(idx));
            D0s = -obj.zb_s(idx);
            D0n = -obj.zb_n(idx);

            nb = obj.neighbors(:,idx);
            dom = obj.domains(idx);

            if dom==0 || dom==1 || dom == 5
                dsig = sig_center(nb(2))-sig_center(nb(4)); % central difference
            elseif dom ==2 || dom==3|| dom==4
                dsig = sig_center(idx)-sig_center(nb(4)); % one-sided difference (deprecated)
            else
                dsig = sig_center(nb(2))-sig_center(idx); % one-sided difference (deprecated)
            end

            if dom==0 || dom==3 || dom == 7
                dn = n_center(nb(1))-n_center(nb(3)); % central difference
            elseif dom ==1 || dom==2|| dom==8
                dn = n_center(idx)-n_center(nb(3)); % one-sided difference (deprecated)
            else
                dn = n_center(nb(1))-n_center(idx); % one-sided difference (deprecated)
            end
            if dom == 0
                if strcmp(tid_name, 'M0A')
                    col = [find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',idx,': d^1u/dx^1: M0A'))) , ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',nb(2),': u0: M0A'))) , ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',nb(4),': u0: M0A'))) , ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',nb(1),': v0: M0A'))) , ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',nb(3),': v0: M0A'))) , ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',nb(2),': v0: M0A'))) , ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',nb(4),': v0: M0A'))) , ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',nb(2),': w0: M0A'))) , ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',nb(4),': w0: M0A')))];
                    row = zeros(size(col));
                    val = [nscale*D0, nscale*(1-sig_center(idx))*D0s/dsig, -nscale*(1-sig_center(idx))*D0s/dsig,...
                        sscale*D0/dn, -sscale*D0/dn, sscale*(1-sig_center(idx))*D0n/dsig, -sscale*(1-sig_center(idx))*D0n/dsig,...
                        sscale*nscale/dsig, -sscale*nscale/dsig];
                else % Tidal constituent introduces two extra terms in the continuity eq.
                    col = [find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',idx,': d^1u/dx^1: ',tid_name))) , ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',nb(2),': u0: ',tid_name))) , ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',nb(4),': u0: ',tid_name))) , ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',nb(1),': v0: ',tid_name))) , ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',nb(3),': v0: ',tid_name))) , ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',nb(2),': v0: ',tid_name))) , ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',nb(4),': v0: ',tid_name))) , ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',nb(2),': w0: ',tid_name))) , ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s%s','cell ',nb(4),': w0: ',tid_name))), ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',idx,': d^1u/dx^1: M0A'))), ...
                        find(strcmp(par_names_tot, sprintf('%s%i%s','cell ',idx,': d^1v/dy^1: M0A')))]; % Two extra terms from subtidal eq
                    row = zeros(size(col));
                    val = [nscale*D0, nscale*(1-sig_center(idx))*D0s/dsig, -nscale*(1-sig_center(idx))*D0s/dsig,...
                        sscale*D0/dn, -sscale*D0/dn, sscale*(1-sig_center(idx))*D0n/dsig, -sscale*(1-sig_center(idx))*D0n/dsig,...
                        sscale*nscale/dsig, -sscale*nscale/dsig, nscale*wl.parameters(tid_idx), sscale*wl.parameters(tid_idx)];
                end
            else
                col = []; row = []; val = [];
            end
        end
    end
end


%                 mesh = mesh;
%             wl = water_level;
%
%             rows = []; cols = []; vals = [];
%             row_idx = 0;
%
%
%             for idx = 1:mesh.ncells %rows
%                 for tid_idx = 1:length(wl.names)
%                     tid_name = wl.names{tid_idx}(end-2:end); % Last three letters
%                     [row, col, val] = obj.idx2element_mat(idx, tid_name, solver);
%                     rows = [rows row_idx+row];
%                     cols = [cols col];
%                     vals = [vals val];
%                     row_idx = max(rows) + 1;
%                 end
%             end