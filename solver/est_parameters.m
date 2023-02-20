function dat = est_parameters(dat, est_opts)

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