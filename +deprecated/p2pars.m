function pars = p2pars(p, ncells)

[np, ne] = size(p);
npars = np/ncells;
pars = zeros(ncells, npars, ne); % could be done using one reshape() call
for n = 1:ne
    pars(:,:,n) = reshape(p(:,n) ,[npars, ncells])';
end
end