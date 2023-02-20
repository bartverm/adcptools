function nullSpace = generate_null_intersection(mat_cell)

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