classdef KinematicRegularization < TaylorBasedRegularization
    methods(Access = protected)
        function assemble_matrix_private(obj)
            assemble_matrix_private@TaylorBasedRegularization(obj);
            if ~obj.model_is_taylor
                return
            end
            if ~(...
                    all(obj.model.sigma_order > 0) ||...
                    all(obj.model.z_order > 0) ...
                )
                warning('KinematicRegularization:TaylorOrderTooLow',...
                    ['No kinematic boundary condition matrix ',...
                    'assembled: Include higher order Taylor expansion',...
                    'in sigma/z direction'])
                return
            end

            % Function that assembles cell-based kinematic boundary conditions.

            const_names = obj.get_const_names(); % Cell array
            zbs = obj.zbsn(1,:);
            zbn = obj.zbsn(2,:);
            dom = obj.domains;
            col = cell([1, numel(const_names)]);
            pnames = obj.flatten_names();
            for eq = 1:numel(const_names)
                col{eq} = [find(strcmp(pnames, ['u0', const_names{eq}])) , ...
                    find(strcmp(pnames, ['d^1u/dsig^1', const_names{eq}])) , ...
                    find(strcmp(pnames, ['v0', const_names{eq}])) , ...
                    find(strcmp(pnames, ['d^1v/dsig^1', const_names{eq}])) , ...
                    find(strcmp(pnames, ['w0', const_names{eq}])) , ...
                    find(strcmp(pnames, ['d^1w/dsig^1', const_names{eq}]))];
            end
            % Now, the indices of the relevant terms are known. Proceed to fill in
            Cj = cell([obj.mesh.ncells,1]);
            rhsj = cell([obj.mesh.ncells,1]);

            for cell_idx = 1:obj.mesh.ncells
                Cj{cell_idx} = zeros([numel(const_names),sum(obj.model.npars)]);
                rhsj{cell_idx} = zeros([numel(const_names),1]);

                if dom(cell_idx)==0 || dom(cell_idx)==1 || dom(cell_idx) == 5 % Internal cells (in terms of sigma)
                    % Do nothing -> No lateral boundary conditions
                else
                    bot = (dom(cell_idx) ==6 || dom(cell_idx)==7|| dom(cell_idx)==8);  %Bottom cells
                    surf = (dom(cell_idx) ==2 || dom(cell_idx)==3|| dom(cell_idx)==4); %Surface cells
                    if surf
                        dsig = 1 - obj.mesh.sig_center(cell_idx); % per definition of sigma. Positive
                        terms = [0, 0, 0, 0, 1, dsig];
                        rhsj{cell_idx} = obj.water_level2rhs_element;
                    elseif bot % Bottom cells
                        dsig = obj.mesh.sig_center(cell_idx); %
                        terms = [zbs(cell_idx), -dsig*zbs(cell_idx), zbn(cell_idx), -dsig*zbn(cell_idx), -1, dsig];
                    end
                    for eq = 1:numel(const_names)
                        Cj{cell_idx}(eq, col{eq}) = terms;
                    end
                end
            end
            obj.C = helpers.spblkdiag(Cj{:});
            obj.rhs = sparse(cell2mat(rhsj));
        end
    end
end