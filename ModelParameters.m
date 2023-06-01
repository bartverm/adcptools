classdef ModelParameters < handle
    %MODELPARAMETERS Class for capturing model parameters as output of
    %get_parameters

    %   Detailed explanation goes here

    properties
        M (:,:) double % matrix M such that b = Mp

        b (:,1) double; %rhs of system of eqs (b = Mp)

        p (:,:) double; % model parameters (b = Mp)

        pars double;

        cov_pars (:,:,:) double;

        regularization (1,1) Regularization

        opts (1,1) SolverOptions

        ns (:,:) double

        cell_idx (:,1) cell

        vel_cmap = brewermap(20, 'RdBu');

        s_cmap = brewermap(15, 'YlOrBr');

        cv_results (1,:) double
    end
    properties(Dependent)
        phi_cmap
    end
    methods
        function obj = ModelParameters(varargin)
            % Overwrite default options
            for ia = 1:2:nargin
                obj.(varargin{ia}) = varargin{ia+1};
            end
        end

        function plot_solution(obj, names_selection, par_idx, varargin)

            if nargin < 2
                names_selection = [obj.regularization.model.names{:}];
            end
            inp = inputParser;
            %inp.addOptional('var',[]);
            inp.addParameter('v',0,@(x) isscalar(x) && isfinite(x));
            inp.addParameter('w',0,@(x) isscalar(x) && isfinite(x));
            expectedTransform = {'linear','symlog'};
            inp.addParameter('ArrowScaling',[.1, .1]);
            inp.addParameter('ArrowTransform','linear', @(x) any(validatestring(x,expectedTransform)));
            inp.addParameter('ArrowParam', [.9, .9])
            inp.parse(varargin{:})
            v = inp.Results.v;
            w = inp.Results.w;
            ArrowTransform = inp.Results.ArrowTransform;
            ArrowScaling = inp.Results.ArrowScaling;
            ArrowParam = inp.Results.ArrowParam;
            P = obj.p(:, par_idx);
            nreg = size(P, 2);
            nn = length(names_selection);
            nc = obj.regularization.mesh.ncells;
            np = sum(obj.regularization.model.npars); % Number of parameters in each cell
            Np = size(P,1); %= nc*np;

            if nreg>1 % Compare different vectors
                t = tiledlayout(nn, nreg, TileSpacing = "tight", Padding = "tight", TileIndexing = "columnmajor");
            else
                t = tiledlayout('flow', TileSpacing="tight", Padding="tight");
            end

            t.XLabel.String = 'y [m]';
            t.YLabel.String = 'z [m]';
            t.XLabel.Interpreter = 'latex';
            t.YLabel.Interpreter = 'latex';

            par_idx = obj.get_par_idx(names_selection);
            titles = obj.modify_names(names_selection);

            for col = 1:nreg
                for row = 1:nn
                    nexttile();
                    amax = max(abs(P(par_idx(row):np:Np, 2)), [], 'omitnan') + 1e-5; % For paper, we assume nreg = 3
                    if ~v && ~w
                        var = [P(par_idx(row):np:Np, col), zeros(nc,2)];
                        armax = [0, 0];
                    elseif v && ~w
                        var = [P(par_idx(row):np:Np, col), P(par_idx(row)+np/3:np:Np, col), zeros(nc,1)];
                        armax = [max(abs(P(par_idx(row)+np/3:np:Np, 2)), [], 'omitnan') + 1e-5, 0];
                    elseif w && ~v
                        var = [P(par_idx(row):np:Np, col), zeros(nc,1), P(par_idx(row)+2*np/3:np:Np, col)];
                        armax = [0, max(abs(P(par_idx(row)+2*np/3:np:Np, 2)), [], 'omitnan') + 1e-5];
                    elseif v && w
                        var = [P(par_idx(row):np:Np, col), P(par_idx(row)+np/3:np:Np, col), P(par_idx(row)+2*np/3:np:Np, col)];
                        armax = [max(abs(P(par_idx(row)+np/3:np:Np, 2)), [], 'omitnan') + 1e-5, max(abs(P(par_idx(row)+2*np/3:np:Np, 2)), [], 'omitnan') + 1e-5];
                    end
                    var = obj.arrow_scale(var, ArrowScaling.*armax, ArrowTransform, ArrowParam);
                    %plot(var(:,2:3))
                    %ylim([-armax(1), armax(1)])
                    hold on
                    obj.regularization.mesh.plot(var, 'FixAspectRatio', false)
                    %     loc_tit = str,old,new)
                    lam = {'\lambda_0', '\lambda_1', '\lambda_2'};
                    title(['$', lam{col},': ' titles{row}, '$'], 'interpreter', 'latex', 'FontSize', 12);

                    %         set(c,'TickLabelInterpreter','latex')
                    % colorbar
                    if col == nreg
                        if ~contains(titles{row}, 'phi')
                            %amax = max(abs(var(:,1)), [], 'omitnan') + 1e-5;
                            c = colorbar;
                            if ~contains(titles{row}, '\partial')
                                ylabel(c, '$m/s$','Rotation',270, 'interpreter', 'latex');
                            elseif contains(titles{row}, '\sigma')
                                ylabel(c, '$m/s$','Rotation',270, 'interpreter', 'latex');
                            else
                                ylabel(c, '$m/s^2$','Rotation',270, 'interpreter', 'latex');
                            end

                        else
                            c = colorbar;
                            ylabel(c, 'deg','Rotation',270, 'interpreter', 'latex');

                        end
                        pos = get(c,'Position');
                        if row == 1
                            pos1 = pos;
                        end
                        %         disp(pos)
                        c.Label.Position(1) = pos1(1)+.5/col; % to change its position
                        %         c.Label.Position(2) = c.Label.Position(2) + .2; % to change its position
                        c.Label.HorizontalAlignment = 'center'; % to change its position
                        c.TickLabelInterpreter = 'latex';
                        %         c.Label.Rotation = 270; % to rotate the text
                    end

                    % colormap
                    if ~contains(titles{row}, 'phi')
                        %amax = max(abs(var(:,1)), [], 'omitnan') + 1e-5;
                        caxis([-amax, amax])

                        colormap(gca, obj.vel_cmap)

                    else
                        caxis([-180, 180])
                        temp = get(gca, 'Children');
                        temp(2).CData =  temp(2).CData*180/pi;
                        colormap(gca, obj.phi_cmap)
                    end

                    axis tight
                    %     hAxes.TickLabelInterpreter = 'latex';
                    %     title(sprintf('%s, %s, %s', names{i}, names{i+np(1)}, names{i + np(1) + np(2)}))
%                     set(gca, 'XDir','reverse') % Very important
                    %         xlabel('y [m]', 'interpreter', 'latex')
                    %         ylabel('z [m]', 'interpreter', 'latex')
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                end
            end
            % Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
            % %[Left Bottom Right Top] spacing
            % NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
            % set(gca, 'Position', NewPos);
        end

        function scaled_var = arrow_scale(obj, var, as, at, ap)%, ArrowScaling, ArrowTransform, ArrowParam);
            % Input: N x 3 array, acts on the last two columns
            scaled_var = var;
            if strcmp(at, 'linear')
                scaled_var(:,2:3) = [var(:,2).*as(1), var(:,3).*as(2)];
            elseif strcmp(at, 'symlog')
                scaled_var(:,2:3) = [helpers.symlog(var(:,2), ap(1)).*as(1), helpers.symlog(var(:,3), ap(2)).*as(2)];
            end
        end

        function par_idx = get_par_idx(obj, names_selection)
            par_idx = nan([1,numel(names_selection)]);
            for name_idx = 1:numel(names_selection)
                par_idx(name_idx) = find(strcmp([obj.regularization.model.names{:}], names_selection{name_idx}));
            end
        end

        function mod_names = modify_names(obj, names_selection)
            mod_names = names_selection;
            for idx = 1:numel(names_selection)
                mod_names{idx} = strrep(mod_names{idx}, 'sig', '\sigma');
                mod_names{idx} = strrep(mod_names{idx}, '^1', '');
                mod_names{idx} = strrep(mod_names{idx}, 'u0', 'u_0');
                mod_names{idx} = strrep(mod_names{idx}, 'v0', 'v_0');
                mod_names{idx} = strrep(mod_names{idx}, 'w0', 'w_0');
                mod_names{idx} = strrep(mod_names{idx}, 'd', '\partial ');
            end
        end


        function training_idx = split_dataset(obj)
            ci = vertcat(obj.cell_idx{:});
            training_idx = ones(size(ci));

            if strcmp(obj.opts.cv_mode, 'none')
                training_idx = ones(size(ci));
            elseif strcmp(obj.opts.cv_mode, 'random')
                tp = obj.opts.training_perc;
                rand0 = rand(size(training_idx));
                training_idx = (rand0 <= tp);
            elseif strcmp(obj.opts.cv_mode, 'omit_cells')
                oc = obj.opts.omit_cells;
                for occ = 1:length(oc)
                    training_idx(ci == oc(occ)) = 0;
                end
            elseif strcmp(obj.opts.cv_mode, 'omit_time') % to be implemented
            end
        end

        function CV = cross_validate_single(obj, reg_pars_cell, pguess)
            % First loop: cross-validation ensembles
            
            Cg = obj.regularization.Cg;
            p_train = cell(size(reg_pars_cell)); %zeros([size(obj.p, 1), nepochs]);
            p_avg = cell(size(reg_pars_cell));

            if strcmp(obj.opts.cv_mode, 'random')
                nepochs = obj.opts.cv_iter;
            else
                nepochs = 1;
            end
            niter = numel(reg_pars_cell)*nepochs;
            i = 0;
            fprintf('Cross-validation percentage: %2.2f percent \n', 100*i/niter)
            for ep = 1:nepochs % randomized iterations: Only meaningful if cv_mode == 'random'
                train_idx = logical(obj.split_dataset());
                test_idx = ~train_idx;

                % Construct training matrix and data
                M0 = obj.M(train_idx, :);
                b0 = obj.b(train_idx);

                M1 = obj.M(test_idx,:);
                b1 = obj.b(test_idx);

                Mp = M0'*M0;

                for rp = 1:numel(reg_pars_cell)
                    regp = reg_pars_cell{rp};
                    A = Mp + regp(1)*Cg{1} + regp(2)*Cg{2} + regp(3)*Cg{3} + regp(4)*Cg{4} + regp(5)*Cg{5};
                    obj.opts.preconditioner_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
                    L = ichol(A, obj.opts.preconditioner_opts);
                    [p_train{rp}(:, ep), ~, ~, it] = pcg(A, M0'*b0 + regp(5)*obj.regularization.C{5}'*obj.regularization.rhs, 1e-9, size(A,2), L, L', pguess(:,rp)); % Matrix of solutions (columns) belonging to regularization parameters regP (rows)
                    i = i+1;
                    fprintf('Cross-validation percentage: %2.2f percent \n', 100*i/niter)
                end
            end
            for rp = 1:numel(reg_pars_cell)
                p_avg{rp} = mean(p_train{rp}, 2);
                CV{rp} = mean((M1*p_avg{rp} - b1).^2);
                disp(['Cross-validated generalization error: ', num2str(CV{rp})])
            end
        end



        function CV = cross_validation(obj)
            % First loop: cross-validation ensembles
            CV = obj.cross_validate_single(obj.opts.reg_pars, obj.p);
        end

        %function solve_single()

        function CV = cross_validation_analysis(obj)
            reg_pars_sens_vec = obj.opts.vectorize_reg_pars();
            for rp = 1:size(reg_pars_sens_vec,1)
                reg_pars_cell{rp} = reg_pars_sens_vec(rp,:);
            end
            CV = obj.cross_validate_single(reg_pars_cell, zeros(size(obj.p,1), numel(reg_pars_cell)));
        end


        function pars = p2pars(obj)
            [np, ne] = size(obj.p);
            ncells = obj.regularization.mesh.ncells;
            npars = np/ncells;
            pars = zeros(ncells, npars, ne); % could be done using one reshape() call, %TODO chatGPT
            for n = 1:ne
                pars(:,:,n) = reshape(obj.p(:,n) ,[npars, ncells])';
            end
            obj.pars = pars;
        end

        function p = pars2p(obj)
            p = reshape(obj.pars', 1, [])';
            obj.p = p;
        end
    end
end


        %         function SE = local_sensitivity(obj)
        %         end
        %
        %         function SE = sensitivity_analysis(obj)
        %
        %         end

        %         function phi_cmap = get.phi_cmap(obj)
        %         phi_cmap = [obj.vel_cmap;  flipud(obj.vel_cmap)];
        %         end

        %     end



% Bin
%                 obj.velocity_model.get_parameter_names();
%
%                 % Internal continuity matrix
%                 %                 for j = 1:obj.mesh.ncells
%                 C1 = assembleC1(obj);
%                 %                 end
%                 C1p = C1'*C1;
%
%                 % External continuity matrix
%                 C2 = assembleC2(obj);
%                 C2p = C2'*C2;
%
%                 %Coherence matrix
%                 C3 = assembleC3(obj);
%                 C3p = C3'*C3;
%
%                 %Consistency matrix
%                 C4 = assembleC4(obj);
%                 C4p = C4'*C4;
%
%                 % Kinematic boundary condition matrix
%                 [C5, bc] = assembleC5(obj);
%                 C5p = C5'*C5;
%
%
%                 % Data matrix: All data
%                 [Mu, Mv, Mw] = obj.velocity_model.get_model(...
%                     dt, cur_s, dn, dz, dsig);
%                 Mb0 = [Mu.*cur_xform(:,1), Mv.*cur_xform(:,2), Mw.*cur_xform(:,3)]; %Model matrix times unit vectors q
%
%                 % If cross-validation is to be applied: Split data in two
%                 % sets. All other operations before are data-independent
%                 % and can thus be performed only once.
%                 Mj = cell([obj.mesh.ncells,1]); ns = zeros([obj.mesh.ncells,1]); bj = cell([obj.mesh.ncells,1]);
%                 for j = 1:obj.mesh.ncells
%                     Mj{j} = Mb0(j==cell_idx,:);
%                     ns(j) = sum(j == cell_idx);
%                     bj{j} = cur_vel(j==cell_idx);
%                 end
%
%                 b = vertcat(bj{:});
%                 M = spblkdiag(Mj{:});
%                 pcg_opts = struct('michol','on','type','ict','droptol',1e-3,'diagcomp',0);
%
%                 % Generate first guess for p
%                 rpg0 = [opts.reg_pars0{:}];
%                 A = M'*M + rpg0(1)*C1p + rpg0(2)*C2p + rpg0(3)*C3p + rpg0(4)*C4p + rpg0(5)*C5p;
%                 alpha = max(sum(abs(A),2)./diag(A))-2;
%                 pcg_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
%                 L = ichol(A, pcg_opts);
%                 p0 = pcg(A, M'*b + rpg0(5)*C5'*bc, 1e-9, size(A,2), L, L'); % global, iterative inversion
%                 np = length(p0);
%
%                 % Compute quantities of interest for comparison
%                 rpg1 = opts.reg_pars1;
%                 p1 = nan(np, size(rpg1,1));
%                 for n = 1:length(rpg1{1})
%                     A = M'*M + rpg1{1}(n)*C1p + rpg1{2}(n)*C2p + rpg1{3}(n)*C3p + rpg1{4}(n)*C4p + rpg1{5}(n)*C5p;
%                     alpha = max(sum(abs(A),2)./diag(A))-2;
%                     pcg_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
%                     p1(:,n) = pcg(A, M'*b + rpg1{5}(n)*C5'*bc, 1e-9, size(A,2), L, L'); % global, iterative inversion
%                 end

%                assignin("base", "dat", struct('M', M, 'C1', C1, 'C2', C2, 'C3', C3, 'C4', C4, 'C5', C5, 'bc', bc, ...
%                    'IM', IM, 'p', p, 'p0', p0, 'p1', p1, 'b', b,...
%                   'opts', opts, 'cell_idx', cell_idx, 'Pe', Pe, 'regP', regP))

%               pars{1,1} = reshape(squeeze(p(:,1,1)) ,[size(Mb0,2), obj.mesh.ncells])'; % At this stage, pars are already known.
%                cov_pars{1,1} = 0; n_vels{1,1} = ns;


            %Continue here.
            %
            %             obj.cv_results
            %             for rp = 1:numel(obj.opts.reg_pars)
            %                 obj.cv_results(rp) = mean((b1 - p_avg).^2)
            % Residuals and goodness of fit

            %                     pe(1, rp, ep) = calc_res(b, M*p(:, rp, ep)); % Performance on full set
            %                     pe(2, rp, ep) = calc_res(b0, M0*p(:, rp, ep)); % Performance on training set
            %                     pe(3, rp, ep) = calc_res(b1, M1*p(:, rp, ep)); % Performance on validation set
            %                     pe(4, rp, ep) = calc_res(0, C1*p(:, rp, ep)); % Performance on continuity
            %                     pe(5, rp, ep) = calc_res(0, C2*p(:, rp, ep)); % Performance on gen. continuity
            %                     pe(6, rp, ep) = calc_res(0, C3*p(:, rp, ep)); % Performance on smoothness
            %                     pe(7, rp, ep) = calc_res(0, C4*p(:, rp, ep)); % Performance on consistency
            %                     pe(8, rp, ep) = calc_res(bc, C5*p(:, rp, ep)); % Performance on boundary conditions
            %
            %                     pe(9,rp,ep) = condest(A);
            %                     pe(10,rp,ep) = it;
