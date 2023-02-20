function p = pars2p(pars)

pars = pars';
p = pars(:);%reshape(squeeze(p(:,1,1)) ,[size(Mb0,2), obj.mesh.ncells])';

end