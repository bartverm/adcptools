function RP = combine_regpars(opts)

% reg_pars must be a cell array of column vectors of not necessarily same size

if strcmp(opts.reg_vary, 'coupled')
    L = opts.reg_pars(1, [1,3]); % Only vary two reg. parameters
    n = length(L);
    [L{:}] = ndgrid(L{end:-1:1});
    L = cat(n+1,L{:});
    L = fliplr(reshape(L,[],n));
    RP = [L(:, 1), L(:,1), L(:,2), L(:,2), L(:,1)]; % Coupling of parameters
elseif strcmp(opts.reg_vary, 'decoupled')
    L = opts.reg_pars; % Vary all reg. parameters
    n = length(L);
    [L{:}] = ndgrid(L{end:-1:1});
    L = cat(n+1,L{:});
    RP = fliplr(reshape(L,[],n));
elseif strcmp(opts.reg_vary, 'none')
    RP = [opts.reg_pars0{:}];
end

for i = 1:size(RP,2)
    if opts.set_zero(i)
        RP(:,i) = 0;
    end
end
end