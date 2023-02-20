function res = matdivide1(A, B)
% Matrix division pointwise while accounting for Inf and Nan
res = A./B;
res(res == Inf) = 1e16; %1 / 0 = 1e16
res(isnan(res)) = 1; % 0 / 0 = 1

end