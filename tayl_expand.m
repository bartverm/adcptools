function [mu,mv,mw]=tayl_expand(ordersu,ordersv,ordersw,in)
%TAYL_EXPAND computes taylor expansions for use with procTrans
%
%   [mu,mv,mw]=TAYL_EXPAND(ORDERSU,ORDERSV,ORDERSW,IN) to compute the
%       Taylor expansions of the optional variables. The optional variables
%       should all be column-vectors of equal size. ORDERSU,ORDERSV and 
%       ORDERSW are vectors with the same number of elements as the input 
%       variables. Since procTrans passes 5 variables (n,z,Sig,s,t) these 
%       variable should have 5 elements indicating to what order the Taylor
%       expansion should be done. A -1 indicates no expansion is performed.
% 
%   Example:
%       orderu=[0 0 -1 1 -1]; %To fit u constant in n and z, and first
%                             % order to s
%       orderv=orderu; % Do the same for v
%       orderw=orderu; % Do the same for w
%
%       procTrans(adcp,tid,'Model',...
%                  @(n,z,Sig,s,t)...
%                  tayl_expand(ordersu,ordersv,ordersv,[n,z,Sig,s,t]),...
%                  'GetVelocity',...
%                  @(pars) tayl_predict(ordersu,ordersv,ordersv,pars))


assert(isequal(numel(ordersu), numel(ordersv), numel(ordersw), size(in,2)))
assert(isequal(ordersu,fix(ordersu)) && all(ordersu>=-1))
assert(isequal(ordersv,fix(ordersv)) && all(ordersv>=-1))
assert(isequal(ordersw,fix(ordersw)) && all(ordersw>=-1))

if isempty(in)
    mu=double.empty(0,npars(ordersu));
    mv=double.empty(0,npars(ordersv));
    mw=double.empty(0,npars(ordersw));
    return
end

mu=expand_var(in,ordersu);
mv=expand_var(in,ordersv);
mw=expand_var(in,ordersw);

end

function n=npars(x)
    n= sum(x+1)-sum(x>-1)+double(any(sum(x>-1)));
end

function out=expand_var(var,orders)
np=npars(orders);
nvels=size(var,1);
out=nan(nvels,np);
cc=0;
has_const=0;
for cp=1:numel(orders)
    for ct=has_const:orders(cp)
        has_const=1;
        cc=cc+1;
        out(:,cc)=(1/factorial(ct))*var(:,cp).^ct;
    end
end

end