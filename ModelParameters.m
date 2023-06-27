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

        GOF (1,:) struct % gof: goodness of fit measures.

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

        function ax = plot_mesh(obj, varargin)

            ax = obj.regularization.mesh.plot(varargin{:});

            xlabel('y [m]')
            if any(strcmp(varargin,'sig')) %dirty
                ylabel('$\sigma$', 'Interpreter', 'latex')
            else
                ylabel('$z [m]$', 'Interpreter', 'latex')
            end
            axis tight
            set(ax, 'XDir','reverse')
            ax.TickLabelInterpreter = 'latex';
        end

        function [me, vare, mse, e] = get_residual(obj, A, x, b)
            % Input: matrix A, solution x, data b
            % Ax = b + e
            % OLS estimator x will produce mimimal residual norm
            % Biased estimators will produce larger residual norms
            % Output: [me, vare, mse, e] = [mean(e), sample variance e,
            % mse(e), e]
            % variance computed as 1/(n-1)*sum(e_i - mean(e))^2

            e = A*x-b;
            me = mean(e);
            mse = e'*e./(numel(e)-1);
            if isnan(me)
                mse = nan;
            end
            vare = mse - numel(e)*me^2/(numel(e)-1);
        end


        function idx = eq2mesh(obj, neq)
            % Mapping between equation index and
            % mesh cell index

            % neq(cell_idx) = number of equations within cell_idx
            % idx{cell_idx} = equation indices within cell_idx
            idx = cell([numel(neq),1]);
            cur_idx = 1;
            for cidx = 1:numel(neq)
                idx{cidx} = [cur_idx:(cur_idx + neq(cidx) - 1)]';
                cur_idx = cur_idx + neq(cidx);
            end
        end

        function [me, vare, mse, e] = get_residual_mesh(obj, A, x, b, neq)
            % Function that obtains the residuals for each mesh cell.
            % Provide the large matrix, solution vector, and rhs
            % Also provide the number of eqs per mesh cell

            % Output: same as get_residual, but now splitted out between
            % mesh cells (added dimension of size mesh.ncells)

            % Implementation using for loop, not very efficient

            me = nan([numel(neq),1]);
            vare = nan([numel(neq),1]);
            mse = nan([numel(neq),1]);
            e = cell([numel(neq),1]);

            idx = obj.eq2mesh(neq);

            for cidx = 1:numel(neq)
                [me(cidx,1), vare(cidx,1), mse(cidx,1), e{cidx, 1}] = obj.get_residual(A(idx{cidx}, :), x, b(idx{cidx}, 1));
                disp(me(cidx,1))
            end
        end

        function gof=get_residuals(obj)
            % Data residuals per mesh:
            for n = 1:size(obj.p,2)
                % Data
                [obj.GOF(n).m, obj.GOF(n).var, obj.GOF(n).mse, obj.GOF(n).e] = obj.get_residual(obj.M, obj.p(:,n), obj.b);
                [obj.GOF(n).mm, obj.GOF(n).varm, obj.GOF(n).msem, obj.GOF(n).em] = obj.get_residual_mesh(obj.M, obj.p(:,n), obj.b, obj.ns);
                % Constraints
                for nc = 1:5
                    if nc < 5
                        rhs = zeros([size(obj.regularization.C{nc}, 1),1]);
                    else
                        rhs = obj.regularization.rhs;
                    end
                    [obj.GOF(n).cm{nc}, obj.GOF(n).cvar{nc}, obj.GOF(n).cmse{nc}, obj.GOF(n).ce{nc}]...
                        = obj.get_residual(obj.regularization.C{nc}, obj.p(:,n), rhs);
                    [obj.GOF(n).cmm{nc}, obj.GOF(n).cvarm{nc}, obj.GOF(n).cmsem{nc}, obj.GOF(n).cem{nc}]...
                        = obj.get_residual_mesh(obj.regularization.C{nc}, obj.p(:,n), rhs, obj.regularization.neq{nc});
                end
            end
            gof = obj.GOF;
        end


        function plot_residual(obj, p_idx, meas_name, var_idx)
            % Plots mesh-based residuals with respect do the data and
            % regularization constraints
            % measure_name = subarray of {"m", "var", "mse"}
            % var_idx = 0,1,2,3,4,5 (0: data, 1 - 5: constraints)

            if var_idx == 0
                fi = strcat(meas_name, "m");
                var = obj.GOF(p_idx).(fi);
                %var(var==0) = nan;
                obj.regularization.mesh.plot('var', var, 'FixAspectRatio', false)
                amax = max(abs(var), [], 'omitnan');
                %if strmp(meas_name, m)
                caxis([-amax, amax])
                %end
            else
                fi = strcat("c", meas_name, "m");
                var = obj.GOF(p_idx).(fi){var_idx};
                amax = max(abs(var), [], 'omitnan');
                %var(var==0) = nan;
                obj.regularization.mesh.plot('var', var, 'FixAspectRatio', false)
                caxis([-amax, amax])
            end        
        end

%         function meas_tit = modify_titles(obj, )


        function plot_residuals(obj, var_idx)
            nreg = size(obj.p, 2);
            nn = numel(var_idx);
            
            for fig_idx = 1:3
                m = makefigure(20,3*numel(var_idx));
                if nreg>1 % Compare different vectors
                    t = tiledlayout(nn, nreg, TileSpacing = "tight", Padding = "tight", TileIndexing = "columnmajor");
                else
                    t = tiledlayout('flow', TileSpacing="tight", Padding="tight");
                end

                t.XLabel.String = 'y [m]';
                t.YLabel.String = 'z [m]';
                t.XLabel.Interpreter = 'latex';
                t.YLabel.Interpreter = 'latex';
                var_tit = {'M\vec{p} - \vec{b}', 'C_1\vec{p}',...
                    'C_2\vec{p}', 'C_3\vec{p}', 'C_4\vec{p}', 'C_5\vec{p} - \vec{\zeta}'};
                meas_name = {'m', 'var', 'mse'};
                meas_tit = {'n_s^{-1}\Sigma', 'Var', 'MSE'};
                for col = 1:nreg
                    for row = 1:nn
                        nexttile;
                        obj.plot_residual(col, meas_name{fig_idx}, var_idx{row})

                        lam = {'\mathbf{\lambda}_0', '\mathbf{\lambda}_1', '\mathbf{\lambda}_2'};
                        if row == 1
                            title(['$', meas_tit{fig_idx}, '(',  var_tit{var_idx{row}+1}, ')', ', \hspace{.1cm}  \lambda = ', lam{col}, '$'], 'interpreter', 'latex', 'FontSize', 12);
                        else
                            title(['$', meas_tit{fig_idx}, '(',  var_tit{var_idx{row}+1}, ')', '$'], 'interpreter', 'latex', 'FontSize', 12);
                        end %['$', titles{row}, ', \hspace{.1cm}  \lambda = ', lam{col}, '$']
                        %tit = strcat();
                        %title(tit, 'interpreter', 'latex', 'FontSize', 12);
                        c = colorbar();
                        set(c,'TickLabelInterpreter','latex')
                        colormap(gca, flipud(obj.vel_cmap))
                        c.FontSize = 12;
                        
                        axis tight
                        set(gca, 'XDir','reverse') % Very important

                        set(gca, 'XTick',[])
                        set(gca, 'YTick',[])

                    end
                end
            end

        end

        function plot_solution(obj, names_selection, par_idx, varargin)
            
            if nargin < 2
                names_selection = [obj.regularization.model.names{:}];
            end
            inp = inputParser;
            %inp.addOptional('var',[]);
            inp.addParameter('v', 0, @(x) isscalar(x) && isfinite(x));
            inp.addParameter('w', 0, @(x) isscalar(x) && isfinite(x));
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
            m = makefigure(20, 3*nn);
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
                    obj.regularization.mesh.plot('var', var, 'FixAspectRatio', false)
                    if col == nreg
                        if ~contains(titles{row}, 'phi')
                            %amax = max(abs(var(:,1)), [], 'omitnan') + 1e-5;
                            c = colorbar;
                            set(c,'TickLabelInterpreter','latex')
                        else
                            c = colorbar;
                            ylabel(c, 'deg','Rotation',270, 'interpreter', 'latex');

                        end
                    end
                    lam = {'\mathbf{\lambda}_0', '\mathbf{\lambda}_1', '\mathbf{\lambda}_2'};

                    if ~contains(titles{row}, '\partial') % Velocities
                        unit = '[m/s]';
                    elseif contains(titles{row}, '\sigma') % Velocities
                        unit = '[m/s]';
                    else % derivatives of velocities in x,y,z directions: m/s/m = 1/s
                        unit = '[1/s]';
                    end
                    if row == 1
                        title(['$', titles{row}, ', \hspace{.1cm}  \lambda = ', lam{col}, '$'], 'interpreter', 'latex', 'FontSize', 12);
                    else
                        title(['$', titles{row}, '$'], 'interpreter', 'latex', 'FontSize', 12);
                    end

                    if col == nreg
                        pos = get(c,'Position');
                        if row == 1
                            pos1 = pos;
                        end
                        % disp(pos)
                        c.Label.String = unit;
                        c.Label.Interpreter = 'latex';
                        %c.Label.Position(1) = pos1(1) + 3; % to change its position
                        %c.Label.Position(2) = c.Label.Position(2) + .2; % to change its position
                        c.Label.HorizontalAlignment = 'center'; % to change its position
                        c.TickLabelInterpreter = 'latex';
                        c.Label.Rotation = 270; % to rotate the text
                        c.FontSize = 12;
                    end

                    % colormap
                    if ~contains(titles{row}, 'phi')    % linear variable
                        caxis([-amax, amax])
                        colormap(gca, obj.vel_cmap)
                        %                         if col == nreg
                        %                         if amax < 1e-2 % change colorbar ticklabels and ylabel to remove scientific notation
                        %                             c.Ticks = 100*c.Ticks;
                        %                             c.Label.String = [c.Label.String, "$\times 10^{-2}$"];
                        %                         end
                        %                         end
                    else
                        caxis([-180, 180])              % cyclic variable
                        temp = get(gca, 'Children');
                        temp(2).CData =  temp(2).CData*180/pi;
                        colormap(gca, obj.phi_cmap)
                    end

                    axis tight
                    set(gca, 'XDir','reverse') % Very important

                    set(gca, 'XTick',[])
                    set(gca, 'YTick',[])

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

            p_train = cell(size(reg_pars_cell));

            if strcmp(obj.opts.cv_mode, 'random')
                nepochs = obj.opts.cv_iter;
            else
                nepochs = 1;
            end
            niter = numel(reg_pars_cell)*nepochs;
            E = cell([numel(reg_pars_cell), 2]);
            CV = cell([numel(reg_pars_cell), 2]);
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
                    p_train{rp}(:, ep) = obj.assemble_solve_single(M0, b0, Mp, regp, pguess(:,rp));
                    i = i+1;
                    fprintf('Cross-validation percentage: %2.2f percent \n', 100*i/niter)
                    E{rp,1}(1, ep) = mean((M1*p_train{rp}(:, ep) - b1).^2); % Generalization error
                    E{rp,2}(1, ep) = mean((M0*p_train{rp}(:, ep) - b0).^2); % Training error
                end
            end
            for rp = 1:numel(reg_pars_cell)
                %                 %p_avg{rp} = mean(p_train{rp}, 2); % k-fold cross-validated estimate
                CV{rp, 1} = mean(E{rp,1}); % Ensemble average
                CV{rp, 2} = mean(E{rp,2}); % Ensemble average
            end
        end

        function [A, rhs, L] = assemble_single(obj, M, b, Mp, regp)
            A = Mp + regp(1)*obj.regularization.Cg{1} + regp(2)*obj.regularization.Cg{2} +...
                regp(3)*obj.regularization.Cg{3} + regp(4)*obj.regularization.Cg{4} + regp(5)*obj.regularization.Cg{5};
            obj.opts.preconditioner_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
            L = ichol(A, obj.opts.preconditioner_opts);
            rhs = M'*b + regp(5)*obj.regularization.C{5}'*obj.regularization.rhs;
        end

        function p = assemble_solve_single(obj, M, b, Mp, regp, pguess)
            % Solve system of eqs
            [A, rhs, L] = assemble_single(obj, M, b, Mp, regp);
            p = solve_single(obj, A, rhs, L, pguess); % Matrix of solutions (columns) belonging to regularization parameters regP (rows)
        end

        function p = solve_single(obj, A, rhs, L, pguess)
            % Solve system of eqs
            [p, ~, ~, ~] = pcg(A, rhs, obj.opts.pcg_tol, obj.opts.pcg_iter, L, L', pguess); % Matrix of solutions (columns) belonging to regularization parameters regP (rows)
        end

        function rhs = b2rhs(obj, M, b, regp)
            rhs =  M'*b + regp(5)*obj.regularization.C{5}'*obj.regularization.rhs;
        end

        function CV = cross_validation(obj)
            % First loop: cross-validation ensembles
            CV = obj.cross_validate_single(obj.opts.reg_pars, obj.p);
        end

        %function solve_single()

        function CV = cross_validation_analysis(obj)
            reg_pars_cell = reg_pars2cell(obj);
            CV = obj.cross_validate_single(reg_pars_cell, zeros(size(obj.p,1), numel(reg_pars_cell)));
        end

        function SA = local_sensitivity(obj)
            % Perform sensitivity analysis using iid noise ~N(0, sigma^2),
            % only on the state vectors that are estimated in the
            % get_parameters call of Solver.
        end


        function [Pe, P] = local_sensitivity_analysis(obj)
            %             M = dat.M; C1 = dat.C1; C2 = dat.C2; C3 = dat.C3; C4 = dat.C4; C5 = dat.C5; bc = dat.bc; b = dat.b;
            %             Np = size(M,2); % number of parameters
            %             nc = max(dat.cell_idx); % number of cells
            %             np = Np/nc; %number of parameters per cell
            %             opts = dat.opts;
            %             nepochs = opts.cv_iter; %nepochs now serves the role
            %             Mp = M'*M; C1p = C1'*C1; C2p = C2'*C2; C3p = C3'*C3; C4p = C4'*C4; C5p = C5'*C5;
            %             regP = combine_regpars(opts);
            %             pcg_opts = struct('michol','on','type','ict','droptol',1e-3);
            reg_pars_cell = reg_pars2cell(obj);
            if strcmp(obj.opts.generate, 'nullspace') % deprecated
                %     D0 = C3;
                %     for row = 1:size(D0,1) % rescale smoothing matrix back - is this necessary?
                %         D0(row,:) = D0(row,:)/C3(row,row);
                %     end
                %     tic;
                %                 NS = generate_null_intersection({C1, C2, C3, C4, C5});
                %     to = toc;
                %                 fprintf('Finished calculating intersection of null spaces after %2.2f s \n', to)

                %                 p = NS*randn(size(NS,2), nepochs);
            elseif strcmp(obj.opts.generate, 'local')
                Mp = obj.M'*obj.M;
                p0 = nan([size(obj.p,1),numel(reg_pars_cell)]);
                for rp = 1:numel(reg_pars_cell) % generate solution for all regularization parameters, using all data
                    regp = reg_pars_cell{rp};
                    p0(:, rp) = obj.assemble_solve_single(obj.M, obj.b, Mp, regp, zeros([size(obj.p,1),1]));
                    %p0: Solution state vectors belonging to all
                    %regularization parameters in the sensitivity analysis.
                    fprintf('Ensemble generation percentage: %2.2f percent \n', 100*rp/numel(reg_pars_cell))
                end
                %p0 = mean(p,2); % Selection of weighted parameter vectors across many reg.pars.
            end


            W0 = rand([numel(reg_pars_cell), obj.opts.ens_size]);
            W = W0/diag(sum(W0,1)); %Scale all columns
            P = p0*W; % True vectors: columns of P

            stdn = obj.opts.get_noise();

            B = obj.M*P; % Unperturbed data: ns x nP
            niter = numel(reg_pars_cell)*numel(stdn)*obj.opts.ens_size*obj.opts.sa_iter; % total number of pcg runs
            fprintf('Total number of iterations will be: %i \n', niter)
            i=0;
            Mp = obj.M'*obj.M;
            Pe = cell([numel(reg_pars_cell), numel(stdn)]);
            %E = cell([numel(reg_pars_cell), numel(stdn)]);
            for rp = 1:numel(reg_pars_cell)
                regp = reg_pars_cell{rp};
                [A, ~, L] = obj.assemble_single(obj.M, obj.b, Mp, regp);
                for nn = 1:numel(stdn)
                    for it = 1:obj.opts.sa_iter % Bootstrapping random iterations
                        Bp = B + stdn(nn)*randn(size(B));
                        for ep = 1:obj.opts.ens_size
                            % TODO apply variation in sa_iter: Noise term must change.
                            i = i+1;
                            % Perturbed measurements -> perhaps move to inner loop.
                            rhs = obj.b2rhs(obj.M, Bp(:,ep), regp); % Estimate different base vector per iteration.
                            fprintf('Sensitivity analysis: %2.2f percent \n', 100*i/niter)
                            Pe{rp, nn}(:, ep, it) = obj.solve_single(A, rhs, L, P(:,ep));
                        end
                    end
                end
            end
        end

        function E = extract_sensitivity_data(obj, Pe, P)
            for rp = 1:size(Pe, 1)
                for nn = 1:size(Pe, 2)
                    for it = 1:obj.opts.sa_iter % Bootstrapping random iterations
                        for ep = 1:obj.opts.ens_size
                            E{rp, nn, 1}(ep, it) = mean((P(:, ep) - Pe{rp, nn}(:, ep, it)).^2); %MSE
                            E{rp, nn, 2}(ep, it) = mean((P(:, ep) - Pe{rp, nn}(:, ep, it)).^2);
                        end
                    end
                end
            end
        end
        %err(rp, ep) =
        %rel_err(rp, ep) = err(rp, ep)./sqrt(mean(P(:,ep).^2));

        % Investigate sensitivity wrt all parameters small quantities
        %                         for ip = 1:np
        %                             perr(ip,rp, ep) = sqrt(sum((p(ip:np:end,ep)-phat(ip:np:end, rp, ep) ).^2));
        %                             rel_perr(ip,rp, ep) = perr(ip,rp, ep)./sqrt(sum(p(ip:np:end,ep).^2));
        %                         end

        %                     avg_err(rp, nn) = mean(err(rp, :));
        %                     avg_perr(:, rp, nn) = mean(squeeze(perr(:, rp, :)),2);
        %                     avg_rel_err(rp, nn) = mean(rel_err(rp, :));
        %                     avg_rel_perr(:, rp, nn) = mean(squeeze(rel_perr(:, rp, :)),2);



        function reg_pars_cell = reg_pars2cell(obj)
            reg_pars_sens_vec = obj.opts.vectorize_reg_pars();
            reg_pars_cell = cell(size(reg_pars_sens_vec,1), 1);
            for rp = 1:size(reg_pars_sens_vec,1)
                reg_pars_cell{rp} = reg_pars_sens_vec(rp,:);
            end
        end

        function fig = plot_contourf_fig(obj, var, xvar, yvar, scaled, tit)
            fig = makefigure(24,12);
            t = tiledlayout('flow', TileSpacing="tight", Padding="tight");
            %             t = tiledlayout(2,3, TileSpacing="tight", Padding="tight");
            %s=nexttile;

            t.XLabel.String = '$\lambda_c$';
            t.YLabel.String = '$\lambda_g$';
            t.XLabel.Interpreter = 'latex';
            t.YLabel.Interpreter = 'latex';
            for idx = 1:numel(var)

                ax{idx} = nexttile;
                plot_contourf_ax(obj, ax{idx}, var{idx}, xvar, yvar, scaled{idx}, tit{idx});
            end
        end

        function ax = plot_contourf_ax(obj, ax, var, xvar, yvar, scaled, tit)

            if scaled
                var = helpers.matdivrobust(var, var(1,1));
                caxis([1-max(max(abs(var-1))), 1+max(max(abs(var-1)))]);
            end
            hold on
            contourf(xvar, yvar, var, 50, 'LineStyle','none')
            contour(xvar, yvar, var, [1,1], 'LineColor','k')
            hold off
            %xlabel('lc'); ylabel('lg')
            colormap(gca, flipud(brewermap(50, 'RdBu')));

            title(tit, 'interpreter', 'latex')
            idxt = unique([1, round(size(xvar,2)/4), round(size(xvar,2)/2), round(3*size(xvar,2)/4), round(size(xvar,2))]);
            idyt = unique([1, round(size(yvar,1)/4), round(size(yvar,1)/2), round(3*size(yvar,1)/4), round(size(yvar,1))]);
            ax.XTick = xvar(1,idxt);
            ax.YTick = yvar(idyt,1);
            ax.XTickLabel = round(helpers.symexp(xvar(1,idxt), obj.opts.res_near_zero(1)), 1, 'significant');
            ax.YTickLabel = round(helpers.symexp(yvar(idyt,1), obj.opts.res_near_zero(3)), 1, 'significant');
            ax.TickLabelInterpreter = 'latex';
            c = colorbar();
            c.TickLabelInterpreter = 'latex';
        end


        function plot_cross_validation_analysis(obj, CV, scaled)
            [CV, lc, lg] = prep_sens_plot(obj, CV);

            fig = plot_contourf_fig(obj, CV, lc, lg, {1,1}, {'Generalization error', 'Training error'});


        end

        function [plot_var, lc, lg] = prep_sens_plot(obj, SA)
            plot_var = cell([size(SA,2), 1]);
            for idx = 1:size(SA,2) % number of metrics to be plotted
                plot_var{idx,1} = reshape(cell2mat(SA(:,idx)),...
                    [obj.opts.reg_iter(3), obj.opts.reg_iter(1)]);
                % CV contains metrics belonging to different reg-pars
            end
            if strcmp(obj.opts.reg_vary, 'coupled')
                [lc, lg] = meshgrid(helpers.symlog(obj.opts.reg_pars_sens{1}, obj.opts.res_near_zero(1)),...
                    helpers.symlog(obj.opts.reg_pars_sens{3}, obj.opts.res_near_zero(3)));
            else
                warning("Use SolverOptions.reg_vary = 'coupled' to perform sensitivity analyses.")
            end
            %             plot_var
            %             lc
            %             lg
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

        function nullSpace = generate_null_intersection(obj, mat_cell)
            % Very slow function that takes cell array of (sparse) matrices and computes
            % intersection of their nullspaces, with basis vectors placed in columns of nullSpace
            cur_null = null(full(mat_cell{1}));
            for m = 2:length(mat_cell)
                next_null = null(full(mat_cell{m}));
                pre_null = null([cur_null, -next_null]);
                cur_null = cur_null*pre_null(1:size(cur_null,2),:);
            end
            nullSpace = cur_null;
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
