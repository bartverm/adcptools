classdef KinematicRegularization < Regularization
    properties
        C5
    end
    methods(Access = protected)
        function assemble_matrix_private(obj, bathy, xs, mesh, data_model, water_level)
            assert(isa(data_model,"TaylorModel"), "To impose kinematic BCs, a Taylor data_model is required");
            if (all(data_model.sigma_order > 0) || all(data_model.z_order > 0))
                obj.C5 = assemble_kinematic(obj, solver);
            else
                warning('No kinematic boundary condition matrix assembled: Include higher order Taylor expansion')
            end
        end

        function [C, rhsvec] = assemble_kinematic(obj, bathy, xs, mesh, data_model, water_level)

            % Function that assembles cell-based kinematic boundary conditions.

            mesh = mesh;
            data_model = data_model;
            wl = water_level;
            wl.get_omega();
            % Compute relevant indices

            subtidal_idx = [find(strcmp(data_model.names, sprintf('%s%s%s', 'u0: M0A'))) , ...
                find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1u/dsig^1: M0A'))) ,...
                find(strcmp(data_model.names, sprintf('%s%s%s', 'v0: M0A'))) , ...
                find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1v/dsig^1: M0A'))) , ...
                find(strcmp(data_model.names, sprintf('%s%s%s', 'w0: M0A'))) , ...
                find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1w/dsig^1: M0A')))];
            tidal_idx = cell([numel(data_model.constituentsU),2]);
            
            for i = 1:numel(data_model.constituentsU) %Fundamental choice: All constituents are the same!!
                const = data_model.constituentsU{i};
                tidal_idx{i,1} = [find(strcmp(data_model.names, sprintf('%s%s%s', 'u0: ', const,'A'))) , ...
                    find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1u/dsig^1: ', const,'A'))) ,...
                    find(strcmp(data_model.names, sprintf('%s%s%s', 'v0: ', const,'A'))) , ...
                    find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1v/dsig^1: ', const,'A'))) , ...
                    find(strcmp(data_model.names, sprintf('%s%s%s', 'w0: ', const,'A'))) , ...
                    find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1w/dsig^1: ', const,'A')))];

                tidal_idx{i,2} = [find(strcmp(data_model.names, sprintf('%s%s%s', 'u0: ', const,'B'))) , ...
                    find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1u/dsig^1: ', const,'B'))) ,...
                    find(strcmp(data_model.names, sprintf('%s%s%s', 'v0: ', const,'B'))) , ...
                    find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1v/dsig^1: ', const,'B'))) , ...
                    find(strcmp(data_model.names, sprintf('%s%s%s', 'w0: ', const,'B'))) , ...
                    find(strcmp(data_model.names, sprintf('%s%s%s', 'd^1w/dsig^1: ', const,'B')))];
            end

            % Now, the indices of the relevant terms are known. Proceed to fill in
            el_mat_size  = [1+2*numel(wl.constituents), sum(data_model.npars)];
            Cj = cell([mesh.ncells,1]);
            rhs = cell([mesh.ncells,1]);

            for idx = 1:mesh.ncells
                dom = obj.domains(idx);
                Cj{idx} = zeros(el_mat_size);
                rhs{idx} = zeros([el_mat_size(1),1]);
                
                if dom==0 || dom==1 || dom == 5 % Internal cells (in terms of sigma)
                    % Do nothing -> No lateral boundary conditions
                else
                    bot = (dom ==6 || dom==7|| dom==8);  %Bottom cells
                    surf = (dom ==2 || dom==3|| dom==4); %Surface cells
                    if surf
                        dsig = 1-mesh.sig_center(idx); % per definition of sigma. Positive
                        % ind = zeros([3, numel(names)]);
                        % First: subtidal equation
                        subtidal_tidal_terms = [0, 0, 0, 0, 1, dsig];
                        Cj{idx}(1,subtidal_idx) = subtidal_tidal_terms;
                        rhs{idx}(1,1) = 0;
                        for i = 1:numel(data_model.constituentsU) %Fundamental choice: All constituents are the same!!
                            Cj{idx}(2*i,tidal_idx{i,1}) = subtidal_tidal_terms;
                            Cj{idx}(2*i+1,tidal_idx{i,2}) = subtidal_tidal_terms;
                            rhs{idx}(2*i, 1) = wl.omega(i)*wl.parameters(2*i+1);
                            rhs{idx}(2*i+1, 1) = -wl.omega(i)*wl.parameters(2*i);
                        end
                    elseif bot % Bottom cells
                        dsig = mesh.sig_center(idx); %

                        subtidal_tidal_terms = [obj.zb_s(idx), -dsig*obj.zb_s(idx), obj.zb_n(idx), -dsig*obj.zb_n(idx), -1, dsig];
                        Cj{idx}(1,subtidal_idx) = subtidal_tidal_terms;
                        for i = 1:numel(data_model.constituentsU) %Fundamental choice: All constituents are the same!!
                            Cj{idx}(2*i,tidal_idx{i,1}) = subtidal_tidal_terms;
                            Cj{idx}(2*i+1,tidal_idx{i,2}) = subtidal_tidal_terms;
                        end
                    end
                end
            end
            C = spblkdiag(Cj{:});
            rhsvec = cell2mat(rhs');
        end
    end
end