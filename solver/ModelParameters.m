classdef ModelParameters < handle
    %MODELPARAMETERS Class for capturing model parameters as output of
    %get_parameters

    %   Detailed explanation goes here

    properties
        model (1,1) DataModel = DataModel; %Underlying model for scalar or vector quantities

        M (1,:) double = {}; % matrix M such that b = Mp

        b (:,1) double = []; %rhs of system of eqs (b = Mp)

        p (:,:) double = []; % model parameters (b = Mp)

        C (1,:) cell = {}; % regularization matrices

        cell_idx (:,1) double = [];
    end

    methods
        function obj = ModelParameters(varargin)
            % Overwrite default options
            for ia = 1:2:nargin
                obj.(varargin{ia}) = varargin{ia+1};
            end
        end

        function cross_validation(obj)
            % TODO: FIX



            % If cross-validation is to be applied: Split data in two
            % sets. All other operations before are data-independent
            % and can thus be performed only once.

            if ~strcmp(obj.cv_mode, 'random')
                % Only one partition involved in cv analysis
                train_idx = logical(obj.split_dataset(cell_idx));
            else
                train_idx = nan(length(cell_idx), obj.opts.cv_iter);
                for cvi = 1:obj.opts.cv_iter
                    train_idx(:,cvi) = logical(obj.split_dataset(cell_idx));
                end

            end
            test_idx = ~train_idx;


            for idx = 1:length(obj.opts.reg_pars)
                for tidx = 1:size(train_idx,2)
                    M_train = M(train_idx(:,tidx), :);
                    b_train = b(train_idx(:,tidx));

                    M_test = M(test_idx(:,tidx),:);
                    b_test = b(test_idx(:,tidx));

                    Mp_train = M_train'*M_train;
                    rp = obj.opts.reg_pars{idx};
                    A = Mp_train + rp(1)*C1p + rp(2)*C2p + rp(3)*C3p + rp(4)*C4p + rp(5)*C5p;
                    if obj.opts.set_diagcomp
                        obj.opts.preconditioner_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
                    end
                    L = ichol(A, obj.opts.preconditioner_opts);
                    p(:,idx) = pcg(A, M_train'*b_train + rp(5)*C5'*bc, obj.opts.pcg_tol, obj.opts.pcg_iter, L, L');

                end
                %
                %
                %
                %             % First loop: cross-validation ensembles
                %             opts.cell_idx = cell_idx;
                %
                %             % Regularization parameters
                %
                %
                %             nepochs = opts.cv_iter;
                %             p = zeros([np,size(regP,1),nepochs]);
                %             niter = nepochs*size(regP,1);
                %             i = 0;
                %             if opts.use_p0
                %                 pguess = p0;
                %             else
                %                 pguess = zeros(size(p0));
                %             end
                %
                %
                %
                %
                %
                %             if obj.opts.gen_analysis
                %
                %             pe = zeros([10, size(regP,1), nepochs]);
                %             if ~strcmp(opts.cv_mode, 'none')
                %                 for ep = 1:nepochs
                %                     train_idx = logical(split_dataset(opts));
                %                     test_idx = ~train_idx;
                % regP = combine_regpars(opts);
                %                     % Construct training matrix and data
                %                     M0 = M(train_idx, :);
                %                     b0 = b(train_idx);
                %
                %                     M1 = M(test_idx,:);
                %                     b1 = b(test_idx);
                %
                %                     Mp = M0'*M0;
                %
                %                     for rp = 1:size(regP,1)
                %                         i = i+1;
                %                         fprintf('Cross-validation percentage: %2.2f percent \n', 100*i/niter)
                %                         A = Mp + regP(rp,1)*C1p + regP(rp,2)*C2p + regP(rp,3)*C3p + regP(rp,4)*C4p + regP(rp,5)*C5p;
                %                         pcg_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
                %                         L = ichol(A, pcg_opts);
                %                         [p(:, rp, ep), ~, ~, it] = pcg(A, M0'*b0 + regP(rp,5)*C5'*bc, 1e-9, size(A,2), L, L', pguess); % Matrix of solutions (columns) belonging to regularization parameters regP (rows)
                %
                %                         % Residuals and goodness of fit
                %
                %                         pe(1, rp, ep) = calc_res(b, M*p(:, rp, ep)); % Performance on full set
                %                         pe(2, rp, ep) = calc_res(b0, M0*p(:, rp, ep)); % Performance on training set
                %                         pe(3, rp, ep) = calc_res(b1, M1*p(:, rp, ep)); % Performance on validation set
                %                         pe(4, rp, ep) = calc_res(0, C1*p(:, rp, ep)); % Performance on continuity
                %                         pe(5, rp, ep) = calc_res(0, C2*p(:, rp, ep)); % Performance on gen. continuity
                %                         pe(6, rp, ep) = calc_res(0, C3*p(:, rp, ep)); % Performance on smoothness
                %                         pe(7, rp, ep) = calc_res(0, C4*p(:, rp, ep)); % Performance on consistency
                %                         pe(8, rp, ep) = calc_res(bc, C5*p(:, rp, ep)); % Performance on boundary conditions
                %
                %                         pe(9,rp,ep) = condest(A);
                %                         pe(10,rp,ep) = it;
                %                     end
                %                 end
                %             end
                %             Pe = mean(pe, 3);
                %             assignin("base", "dat", struct('M', M, 'C1', C1, 'C2', C2, 'C3', C3, 'C4', C4, 'C5', C5, 'bc', bc, ...
                %                 'IM', IM, 'p', p, 'p0', p0, 'p1', p1, 'b', b,...
                %                 'opts', opts, 'cell_idx', cell_idx, 'Pe', Pe, 'regP', regP))
                %
                %             pars{1,1} = reshape(squeeze(p(:,1,1)) ,[size(Mb0,2), obj.mesh.ncells])'; % At this stage, pars are already known.
                %             cov_pars{1,1} = 0; n_vels{1,1} = ns;
            end

        end

        function sensitivity_analysis(obj)
            % TODO: FIX
            M = dat.M; C1 = dat.C1; C2 = dat.C2; C3 = dat.C3; C4 = dat.C4; C5 = dat.C5; bc = dat.bc; b = dat.b;
            Np = size(M,2); % number of parameters
            nc = max(dat.cell_idx); % number of cells
            np = Np/nc; %number of parameters per cell
            opts = dat.opts;
            nepochs = opts.cv_iter; %nepochs now serves the role
            Mp = M'*M; C1p = C1'*C1; C2p = C2'*C2; C3p = C3'*C3; C4p = C4'*C4; C5p = C5'*C5;
            regP = combine_regpars(opts);
            pcg_opts = struct('michol','on','type','ict','droptol',1e-3);

            if strcmp(est_opts.generate, 'nullspace')
                %     D0 = C3;
                %     for row = 1:size(D0,1) % rescale smoothing matrix back - is this necessary?
                %         D0(row,:) = D0(row,:)/C3(row,row);
                %     end
                %     tic;
                NS = generate_null_intersection({C1, C2, C3, C4});
                %     to = toc;
                fprintf('Finished calculating intersection of null spaces after %2.2f s \n', to)

                p = NS*randn(size(NS,2), nepochs);
            elseif strcmp(est_opts.generate, 'local')
                p = repmat(dat.p0, 1, nepochs); % Change this to be less biased towards p0 (with corresponding reg pars)
                % Generate ensemble
                for rp = 2:size(regP,1)
                    fprintf('Ensemble generation: %2.2f percent \n', 100*rp/size(regP,1))
                    A = Mp + regP(rp,1)*C1p + regP(rp,2)*C2p + regP(rp,3)*C3p + regP(rp,4)*C4p + regP(rp,5)*C5p;
                    pcg_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
                    L = ichol(A, pcg_opts);
                    p(:, rp) = pcg(A, M'*b + regP(rp,5)*C5'*bc, 1e-6, size(A,2), L, L', dat.p0); % Matrix of solutions (columns) belonging to regularization parameters regP (rows)
                end
                p0 = mean(p,2); % Mean parameter vector across many reg.pars.
            end
            % p = nan([np,size(regP,1),nepochs]);
            P = repmat(p0, 1, nepochs);
            stdn = est_opts.noise_levels;

            B = M*P; % Unperturbed data
            niter = length(stdn)*size(regP,1)*nepochs;
            fprintf('Total number of iterations will be: %i \n', niter)
            i=0;
            avg_err = nan([size(regP,1), length(stdn)]);
            avg_perr = nan([np, size(regP,1), length(stdn)]);
            for nn = 1:length(stdn)
                Bp = B + stdn(nn)*randn(size(B)); % Perturb measurements
                for rp = 1:size(regP,1)
                    A = Mp + regP(rp,1)*C1p + regP(rp,2)*C2p + regP(rp,3)*C3p + regP(rp,4)*C4p + regP(rp,5)*C5p;
                    pcg_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
                    L = ichol(A, pcg_opts);
                    for ep = 1:nepochs
                        i = i+1;
                        fprintf('Sensitivity analysis: %2.2f percent \n', 100*i/niter)
                        phat(:, rp, ep) = pcg(A, M'*Bp(:,ep) + regP(rp,5)*C5'*bc, 1e-6, size(A,2), L, L', p(:,ep)); % Matrix of solutions (columns) belonging to regularization parameters regP (rows)

                        err(rp, ep) = sqrt(mean((P(:,ep)-phat(:, rp, ep) ).^2));
                        rel_err(rp, ep) = err(rp, ep)./sqrt(mean(P(:,ep).^2));

                        % Investigate sensitivity wrt all parameters small quantities
                        for ip = 1:np
                            perr(ip,rp, ep) = sqrt(sum((p(ip:np:end,ep)-phat(ip:np:end, rp, ep) ).^2));
                            rel_perr(ip,rp, ep) = perr(ip,rp, ep)./sqrt(sum(p(ip:np:end,ep).^2));
                        end
                    end
                    avg_err(rp, nn) = mean(err(rp, :));
                    avg_perr(:, rp, nn) = mean(squeeze(perr(:, rp, :)),2);
                    avg_rel_err(rp, nn) = mean(rel_err(rp, :));
                    avg_rel_perr(:, rp, nn) = mean(squeeze(rel_perr(:, rp, :)),2);
                end
            end

            dat.avg_err = avg_err;
            dat.avg_perr = avg_perr;
            dat.avg_rel_err = avg_rel_err;
            dat.avg_rel_perr = avg_rel_perr;
        end



        function training_idx = split_dataset(obj, cell_idx)

            training_idx = ones(size(cell_idx));

            if strcmp(obj.opts.cv_mode, 'none')
                training_idx = ones(size(cell_idx));
            elseif strcmp(obj.opts.cv_mode, 'random')
                tp = obj.opts.training_perc;
                rand0 = rand(size(training_idx));
                training_idx = (rand0 <= tp);
            elseif strcmp(obj.opts.cv_mode, 'omit_cells')
                oc = obj.opts.omit_cells;
                for occ = 1:length(oc)
                    training_idx(cell_idx==oc(occ)) = 0;
                end
            elseif strcmp(obj.opts.cv_mode, 'omit_time') % to be implemented
            end
        end
    end

end