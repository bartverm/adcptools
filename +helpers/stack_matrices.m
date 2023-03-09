function C = stack_matrices(mat_cell)
if ismatrix(mat_cell)
    C = mat_cell;
elseif iscell(mat_cell)
    C = zeros(size(mat_cell{1}));
    for i = 1:length(mat_cell)
        C(:,:,i) = mat_cell{i};
    end
end
end
