function pars = p2pars(p, ncells)


np = length(p)/ncells;
pars = reshape(p ,[np, ncells])';



end