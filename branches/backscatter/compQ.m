%    Copyright 2014 Bart Vermeulen
%
%    This file is part of ADCPTools.
%
%    ADCPTools is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    ADCPTools is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with ADCPTools.  If not, see <http://www.gnu.org/licenses/>.
%
%    input
%    msh       : adcpttools mesh
%    trancsts : trancsts to process (default: all)
% 
%
% TODO zero value boundary conditions should be given at bottom and
%      side walls in combination with linear interpolation,
%      while interpolating the top still with the nearest neighbour method    
%      better: quadratic interpolation with enforced bank angle of say 30 deg at surface
% TODO linear extrapolation of bank position may yield implausible large values
%      if the bottom slope is low or if the GPS coordinates contain an outlier
function [msh Q A] = compQ(msh, trancsts, nfrac, vfield)

    % automatically select all trancsts
    if (nargin() < 2 || isempty(trancsts))
            trancsts = 1:length(msh);
    end

    % relative size of cross section that is used to determination the bank position
    %nfrac=[0.1 0.3 0.3 0.5 0.3 0.1 0.1] ;
    if (nargin()<3 || isempty(nfrac))
            nfrac = 0.2*ones(size(trancsts));
    end

    if (nargin < 4)
	    vfield = 'vel';
    end
    progvelfield = ['prog' vfield];
    sfield = 'cs';
  
    % estimate discharge
    for ct=trancsts(:)'
	[A(ct) p  ] = integrate_area(msh(ct).p, msh(ct).N, msh(ct).Z, trancsts,nfrac(ct));
	[Q(ct) vel] = integrate_discharge(A(ct), p, msh(ct).N, msh(ct).Z, msh(ct).(sfield).(vfield));
	msh(ct).p    = p;
	msh(ct).cs.Q = Q(ct);
	msh(ct).cs.A = A(ct);
    end % for ct

    % compute discharge for progressive
    if (isfield(msh(1),'progvel') && ~isempty(msh(1).progvel))
	for ct=1:length(msh)
		secrot=[ cos(msh(ct).(sfield).dir) sin(msh(ct).(sfield).dir) 0;...
        	        -sin(msh(ct).(sfield).dir) cos(msh(ct).(sfield).dir) 0;...
                        0 0 1];
	for idx=1:size(msh(ct).progvel,3)
		vel = squeeze(msh(ct).(progvelfield)(:,:,idx,:));
		% rotate
		vel1 = vel(:,:,1);
		vel2 = vel(:,:,2);
		vel12 = [vel1(:), vel2(:)]*secrot(1:2,1:2)';
		vel(:,:,1:2) = reshape(vel12,[size(vel,1),size(vel,2),2]);
		Q_ = integrate_discharge(A(ct), msh(ct).p, msh(ct).N, msh(ct).Z, vel);
		msh(ct).progQ(idx) = Q_;
	end % for idx
	end % for ct
    end % if progflag
end % function compQ

% p   : msh mesh structure
% vel : velocity associated with the cells
function [A p] = integrate_area(p, N, Z, transcsts,nfrac)
    % estimate surface area of cells
    p.A=0.5*(sum(bsxfun(@times,p.Z([6 1 5 2],:),[1 -1 1 -1]'),1).*...
                    sum(bsxfun(@times,p.N([2 1],:),[1 -1]'),1)+...
                    sum(bsxfun(@times,p.Z([5 2 4 3],:),[1 -1 1 -1]'),1).*...
                    sum(bsxfun(@times,p.N([3 2],:),[1 -1]'),1));

    % Make cells and determine area for bottom cells
    p.Nbot=[p.nbed(1:2:end-2);
                    p.nbed(2:2:end-1);
                    p.nbed(3:2:end);
                    p.nbed(3:2:end);
                    p.nbed(2:2:end-1);
                    p.nbed(1:2:end-2);
                    p.nbed(1:2:end-2)];
    % TODO no magic numbers
    p.Zbot=        [p.zbed(1:2:end-2);
                    p.zbed(2:2:end-1);
                    p.zbed(3:2:end);
                    p.zbed(3:2:end)*0.96;
                    p.zbed(2:2:end-1)*0.96;
                    p.zbed(1:2:end-2)*0.96;
                    p.zbed(1:2:end-2)];
    
    % put additional virtual cell between last measured cell and bottom 
    p.Abot=0.5*(sum(bsxfun(@times,p.Zbot([6 1 5 2],:),[1 -1 1 -1]'),1).*...
                       sum(bsxfun(@times,p.Nbot([2 1],:),[1 -1]'),1)+...
                       sum(bsxfun(@times,p.Zbot([5 2 4 3],:),[1 -1 1 -1]'),1).*...
                       sum(bsxfun(@times,p.Nbot([3 2],:),[1 -1]'),1));
    
    % Make cells and determine area for bottom cells
    p.Ntop = p.Nbot([1 3 4 6 7],:);
    maxz   = nanmax(p.Z(:));
    p.Ztop = repmat([maxz maxz 0 0 maxz]',[1 size(p.Ntop,2)]);
    p.Atop = (sum(bsxfun(@times,p.Ztop([4 1],:),[1 -1]'),1) ...
              .*sum(bsxfun(@times,p.Ntop([2 1],:),[1 -1]'),1));

    % determine n coordinate of bank by linear extrapolation
    ntr   =round(diff(p.nbed([1 end]))*nfrac/nanmedian(diff(p.nbed)));
    poly  =polyfit(p.zbed(1:ntr),p.nbed(1:ntr),1);
    nleft =poly(2);
%    if ct==5, nleft=-10; end
    poly  =polyfit(p.zbed(end-ntr+1:end),p.nbed(end-ntr+1:end),1);
    nright=poly(2);

    % Make cells and determine area for bottom cells
    p.Nleft = [nleft; p.nbed([1 1])'; nleft];
    p.Zleft = [0 p.zbed(1) 0 0]';
    p.Aleft = -0.5*diff(p.Nleft([1 2]))*p.zbed(1);


    p.Nright = [p.nbed(end) nright p.nbed([end end])]';
    p.Zright = [p.zbed(end) 0 0 p.zbed(end)]';
    p.Aright = 0.5*diff(p.Nright([1 2]))*diff(p.Zright([1 2]));

    % cells at the boundary with unknown depth are NaN
    % their area is set to zero (dry)
    p.A(isnan(p.A)) = 0;
    fdx = isnan(p.Abot);
    p.Abot(fdx)     = 0;
%    velBot(fdx)     = 0;
    fdx = isnan(p.Atop);
    p.Atop(fdx)     = 0;
%    velTop(fdx) = 0;
    p.Aleft(isnan(p.Aleft)) = 0;
    p.Aright(isnan(p.Aright)) = 0;

    % compute area
    A.mid   = sum(p.A);
    A.top   = sum(p.Atop);
    A.bot   = sum(p.Abot);
    A.left  = p.Aleft;
    A.right = p.Aright;
    A.total = A.mid+A.top+ A.bot+A.left + A.right;
end % integrate_area

function [Q vel] = integrate_discharge(A,p,N,Z,vel)
    % ensure backward compatibility to old matlab versions
    % TODO use matlab ver function
    verstr = version();
    if (str2num([verstr(1) verstr(3)]) < 81)
            ifunc = @(x,y) TriScatteredInterp(x,y,'nearest');
    else
            ifunc = @(x,y) scatteredInterpolant(x,y,'nearest','nearest');        
    end

    % Extrapolate velocity and fill in missing velocity
    fgood  = isfinite(vel(:,:,1));
    fbad   = ~fgood;
    fgood3 = cat(3,fgood,false(size(fgood)),false(size(fgood)));
    fbad3  = cat(3,fbad,false(size(fgood)),false(size(fgood)));   
    X      = [N(fgood),Z(fgood)];
        % the triangular interpolator does not extrapolate,
        % therefore a covex hull of zeros is created
        %            TU     = feval(ifunc, [X; 1./sqrt(eps)*[-1, -1; -1, +1; +1 -1; +1, +1]], ...
        %                                  [msh(ct).cs.vel(fgood3); [0 0 0 0]']);

    Nbot = p.Nbot(2,:);
    Zbot = nanmean(p.Zbot([2 5],:),1);
    Ntop = Nbot; % TODO : why not top ?
    Ztop = nanmean(p.Ztop([2 3],:),1);
    Nleft    = nanmean(p.Nleft(1:3,:),1);
    Zleft    = nanmean(p.Zleft(1:3,:),1);
    Nright   = nanmean(p.Nright(1:3,:),1);
    Zright   = nanmean(p.Zright(1:3,:),1);

    % extrapolate top, bottom and side flow
    TU         = feval(ifunc, X, vel(fgood3));
    velBot     = TU([Nbot;   Zbot]');
    velTop     = TU([Ntop;   Ztop]');
    velLeft    = TU([Nleft;  Zleft]');
    velRight   = TU([Nright; Zright]');

    % interpolate invalid cells
    vel(fbad3) = TU([N(fbad) Z(fbad)]);

    % compute discharge
    Q.mid   = sum(vel(p.fgood_3(:,1)).*p.A');
    Q.top   = sum(velTop.*p.Atop');
    Q.bot   = sum(velBot.*p.Abot');
    Q.left  = velLeft.*A.left;
    Q.right = velRight.*A.right;
    Q.total = Q.mid + Q.top + Q.bot + Q.left + Q.right;
end % integrate_discharge

