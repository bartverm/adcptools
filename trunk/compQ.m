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
%    transects : transects to process (default: all)
% 
%
% TODO zero value boundary conditions should be given at bottom and
%      side walls in combination with linear interpolation,
%      while interpolating the top still with the nearest neighbour method    
%      better: quadratic interpolation with enforced bank angle of say 30 deg at surface
% TODO linear extrapolation of bank position may yield implausible large values
%      if the bottom slope is low or if the GPS coordinates contain an outlier
function [msh Q A] = compQ(msh,transects,nfrac)
    vfield = 'vele';
    
    % ensure backward compatibility to old matlab versions
    verstr = version();
    if (str2num([verstr(1) verstr(3)]) < 81)
            ifunc = @(x,y) TriScatteredInterp(x,y,'nearest');
    else
            ifunc = @(x,y) scatteredInterpolant(x,y,'nearest','nearest');        
    end

    % automatically select all transects
    if (nargin() < 2)
            transects = 1:length(msh);
    end

    % relative size of section that is used to determination the bank position
    %nfrac=[0.1 0.3 0.3 0.5 0.3 0.1 0.1] ;
    if (nargin()<3)
            nfrac = 0.2*ones(size(transects));
    end
    
    %% estimate discharge
    for ct=transects(:)'
        % estimate surface area of cells
        msh(ct).p.A=.5*(sum(bsxfun(@times,msh(ct).p.Z([6 1 5 2],:),[1 -1 1 -1]'),1).*...
                        sum(bsxfun(@times,msh(ct).p.N([2 1],:),[1 -1]'),1)+...
                        sum(bsxfun(@times,msh(ct).p.Z([5 2 4 3],:),[1 -1 1 -1]'),1).*...
                        sum(bsxfun(@times,msh(ct).p.N([3 2],:),[1 -1]'),1));

        % Make cells and determine area for bottom cells
        msh(ct).p.Nbot=[msh(ct).p.nbed(1:2:end-2);
                        msh(ct).p.nbed(2:2:end-1);
                        msh(ct).p.nbed(3:2:end);
                        msh(ct).p.nbed(3:2:end);
                        msh(ct).p.nbed(2:2:end-1);
                        msh(ct).p.nbed(1:2:end-2);
                        msh(ct).p.nbed(1:2:end-2)];
        msh(ct).p.Zbot=[msh(ct).p.zbed(1:2:end-2);
                        msh(ct).p.zbed(2:2:end-1);
                        msh(ct).p.zbed(3:2:end);
                        msh(ct).p.zbed(3:2:end)*0.96;
                        msh(ct).p.zbed(2:2:end-1)*0.96;
                        msh(ct).p.zbed(1:2:end-2)*0.96;
                        msh(ct).p.zbed(1:2:end-2)];
        msh(ct).Nbot=msh(ct).p.Nbot(2,:);
        msh(ct).Zbot=nanmean(msh(ct).p.Zbot([2 5],:),1);
        
        msh(ct).p.Abot=.5*(sum(bsxfun(@times,msh(ct).p.Zbot([6 1 5 2],:),[1 -1 1 -1]'),1).*...
                           sum(bsxfun(@times,msh(ct).p.Nbot([2 1],:),[1 -1]'),1)+...
                           sum(bsxfun(@times,msh(ct).p.Zbot([5 2 4 3],:),[1 -1 1 -1]'),1).*...
                           sum(bsxfun(@times,msh(ct).p.Nbot([3 2],:),[1 -1]'),1));
        
        % Make cells and determine area for bottom cells
        msh(ct).p.Ntop=msh(ct).p.Nbot([1 3 4 6 7],:);
        maxz          =nanmax(msh(ct).p.Z(:));
        msh(ct).p.Ztop=repmat([maxz maxz 0 0 maxz]',[1 size(msh(ct).p.Ntop,2)]);
        msh(ct).p.Atop=   (sum(bsxfun(@times,msh(ct).p.Ztop([4 1],:),[1 -1]'),1).*...
                           sum(bsxfun(@times,msh(ct).p.Ntop([2 1],:),[1 -1]'),1));
        msh(ct).Ntop=msh(ct).Nbot;
        msh(ct).Ztop=nanmean(msh(ct).p.Ztop([2 3],:),1);


        % determine n coordinate of bank by linear extrapolation
        ntr   =round(diff(msh(ct).p.nbed([1 end]))*nfrac(ct)/nanmedian(diff(msh(ct).p.nbed)));
        p     =polyfit(msh(ct).p.zbed(1:ntr),msh(ct).p.nbed(1:ntr),1);
        nleft =p(2);
%        if ct==5, nleft=-10; end
        p     =polyfit(msh(ct).p.zbed(end-ntr+1:end),msh(ct).p.nbed(end-ntr+1:end),1);
        nright=p(2);

        % Make cells and determine area for bottom cells
        msh(ct).p.Nleft =[nleft; msh(ct).p.nbed([1 1])'; nleft];
        msh(ct).p.Zleft =[0 msh(ct).p.zbed(1) 0 0]';
        msh(ct).p.Aleft =-0.5*diff(msh(ct).p.Nleft([1 2]))*msh(ct).p.zbed(1);
        msh(ct).Nleft   =nanmean(msh(ct).p.Nleft(1:3,:),1);
        msh(ct).Zleft   =nanmean(msh(ct).p.Zleft(1:3,:),1);

        msh(ct).p.Nright=[msh(ct).p.nbed(end) nright msh(ct).p.nbed([end end])]';
        msh(ct).p.Zright=[msh(ct).p.zbed(end) 0 0 msh(ct).p.zbed(end)]';
        msh(ct).p.Aright=0.5*diff(msh(ct).p.Nright([1 2]))*diff(msh(ct).p.Zright([1 2]));
        msh(ct).Nright  =nanmean(msh(ct).p.Nright(1:3,:),1);
        msh(ct).Zright  =nanmean(msh(ct).p.Zright(1:3,:),1);
        
        % Extrapolate velocity and fill in missing velocity
        fgood  = isfinite(msh(ct).sec.(vfield)(:,:,1));
        fbad   = ~fgood;
        fgood3 = cat(3,fgood,false(size(fgood)),false(size(fgood)));
        fbad3  = cat(3,fbad,false(size(fgood)),false(size(fgood)));   
        X      = [msh(ct).N(fgood),msh(ct).Z(fgood)];
            % the triangular interpolator does not extrapolate,
            % therefore a covex hull of zeros is created
            %            TU     = feval(ifunc, [X; 1./sqrt(eps)*[-1, -1; -1, +1; +1 -1; +1, +1]], ...
            %                                  [msh(ct).sec.vel(fgood3); [0 0 0 0]']);

        % extrapolate top, bottom and side flow
        TU     = feval(ifunc, X, msh(ct).sec.(vfield)(fgood3));
        msh(ct).sec.velBot     = TU([msh(ct).Nbot;msh(ct).Zbot]');
        msh(ct).sec.velTop     = TU([msh(ct).Ntop;msh(ct).Ztop]');
        msh(ct).sec.velLeft    = TU([msh(ct).Nleft;msh(ct).Zleft]');
        msh(ct).sec.velRight   = TU([msh(ct).Nright;msh(ct).Zright]');
        % interpolate invalid cells
        msh(ct).sec.(vfield)(fbad3) = TU([msh(ct).N(fbad) msh(ct).Z(fbad)]);

        % cells at the boundary with unknown depth are NaN
        % their area is set to zero (dry)
        msh(ct).p.A(isnan(msh(ct).p.A)) = 0;
        fdx = isnan(msh(ct).p.Abot);
        msh(ct).p.Abot(fdx)     = 0;
        msh(ct).sec.velBot(fdx) = 0;
        fdx = isnan(msh(ct).p.Atop);
        msh(ct).p.Atop(fdx)     = 0;
%        msh(ct).sec.velTop(fdx) = 0;
        msh(ct).p.Aleft(isnan(msh(ct).p.Aleft)) = 0;
        msh(ct).p.Aright(isnan(msh(ct).p.Aright)) = 0;

        % compute discharge
        msh(ct).sec.Q.mid   = sum(msh(ct).sec.(vfield)(msh(ct).p.fgood_3(:,1)).*msh(ct).p.A');
        msh(ct).sec.Q.top   = sum(msh(ct).sec.velTop.*msh(ct).p.Atop');
        msh(ct).sec.Q.bot   = sum(msh(ct).sec.velBot.*msh(ct).p.Abot');
        msh(ct).sec.Q.left  =  msh(ct).sec.velLeft.*msh(ct).p.Aleft;
        msh(ct).sec.Q.right =  msh(ct).sec.velRight.*msh(ct).p.Aright;
        msh(ct).sec.Q.total =  msh(ct).sec.Q.mid+msh(ct).sec.Q.top ...
                             + msh(ct).sec.Q.bot+msh(ct).sec.Q.left ...
                             + msh(ct).sec.Q.right;

        % compute area
        msh(ct).sec.A.mid   = sum(msh(ct).p.A);
        msh(ct).sec.A.top   = sum(msh(ct).p.Atop);
        msh(ct).sec.A.bot   = sum(msh(ct).p.Abot);
        msh(ct).sec.A.left  =  msh(ct).p.Aleft;
        msh(ct).sec.A.right =  msh(ct).p.Aright;
        msh(ct).sec.A.total =  msh(ct).sec.A.mid+msh(ct).sec.A.top ...
                             + msh(ct).sec.A.bot+msh(ct).sec.A.left ...
                             + msh(ct).sec.A.right;

        Q(ct)=msh(ct).sec.Q;
        A(ct)=mes(ct).sec.A;
    end % for ct
end % function compQ

