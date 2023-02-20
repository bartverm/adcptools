function [nc, sigc] = get_cell_centers(n, sig, m)

% Given a list of n and sigma values (not the same size!!), give two
% modified lists of cell centers according to mesh

idx = m.index(0, .5);

col2cell = m.col_to_cell; % 

nc = m.n_middle(m.col_to_cell);



nidx = m;
nc = NaN;
sigc = NaN;
end