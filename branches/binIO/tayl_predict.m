function varargout=tayl_predict(ordersu, ordersv, ordersw, in)
%TAYL_PREDICT computes taylor expansions for use with procTrans
%
%   [u,v,w]=TAYL_EXPAND(ORDERSU,ORDERSV,ORDERSW,PARS) computes the velocity
%       given the estimated parameters from procTrans. ORDERSU,ORDERSV and 
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



fnu=find(ordersu>-1,1,'first');
fnv=find(ordersv>-1,1,'first');
fnw=find(ordersw>-1,1,'first');
npars=@(x) sum(x+1)-sum(x>-1)+double(any(sum(x>-1)));

nu=npars(ordersu);
nv=npars(ordersv);
if nargout > 0, varargout{1}=in(:,fnu); end
if nargout > 1, varargout{2}=in(:,fnv+nu); end
if nargout > 2, varargout{3}=in(:,fnw+nu+nv); end

end