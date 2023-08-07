function C = kron_modified_mat(A, B)
% Places copies of B in th columns of A after multiplying by
% those columns, i.e.
% Can be vectorized but is not the main bottleneck.
assert((size(A,1)==size(B,1) && ismatrix(A) && ismatrix(B)), 'Design matrix sizes should match')

C = zeros([size(A,1), size(A,2)*size(B,2)]);
idx = 1;
for i = 1:size(A,2)
    col = idx:(idx+size(B,2)-1);
    C(:,col) = repmat(A(:,i), [1, size(B,2)]).*B;
    idx = max(col)+1;
end
end