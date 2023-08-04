function C = kron_modified_cell( A, B)
% Does the same as its matrix equivalent, but now by
% concatenating the elements of cells.
assert(((size(A,1)==1) && (size(B,1) == 1)), 'Only row cells allowed')

% Function is not vectorized, can be done but does not
% influence program performance
C = cell([size(A,1), size(A,2)*size(B,2)]);
idx = 1;
for i = 1:size(A, 2)
    for j = 1:size(B, 2)
        C{1,idx} = [A{i}, B{j}];
        idx = idx + 1;
    end
end
end