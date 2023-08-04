function [pAphi, tid_names] = pab2pAphi(p, names)

np = length(names);

pars = p2pars(p,length(p)/np);

[tid_pars,tid_names] = ab2Aphi(pars, names); % Still in pars layout

pAphi = pars2p(tid_pars);

end
