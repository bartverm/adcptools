function res = matdivrobust(A, B)
% Matrix division pointwise while accounting for Inf and Nan
res = A./B;
if any(isinf(res), 'all')
    res(res == Inf) = 1e16; %1 / 0 = 1e16
    warning('Encountered Inf in pointwise matrix division')
end
    
if any(isnan(res), 'all')
    res(isnan(res)) = 1; % 0 / 0 = 1
    warning('Encountered NaN in pointwise matrix division')
end

end