function [Q, A] = compQ(msh, sections, nfrac)
% Computes discharge through cross sections
%
%   Q = compQ(msh) computes the discharge through the mesh defined in msh
%       and obtained with procTrans. Q is an array with the same size as
%       msh, containing the following fields:
%     Q.mid - The discharge through the measured section
%     Q.top - The discharge through the top of the measured section and the
%             water surface
%     Q.bot - The discharge through the bottom of the mesh and the bed
%     Q.left - Discharge in the left edge of the measurements and the banks
%     Q.right - Discharge in the right edge of the measurements and the
%               banks
%     Q.total - Total discharge
%
%   [Q, A] = compQ(msh) also returns the areas of the different parts of
%   the cross section described above
%
%   compQ(msh, sections) allows to specify which cross sections should be
%   processed. Default is all.
%
%   compQ(msh, sections, nfrac) allows to specify which fraction of the bed
%       must be used to find the bank position by linear extrapolation
%   
%   Please note: this function is a very simplistic approach to discharge
%   estimates, especially for the unmeasured parts.
%
%   see also: procTrans, readDeployment

%    Copyright 2019 Bart Vermeulen
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


% TODO zero value boundary conditions should be given at bottom and
%      side walls in combination with linear interpolation,
%      while interpolating the top still with the nearest neighbour method    
%      better: quadratic interpolation with enforced bank angle of say 30 deg at surface
% TODO linear extrapolation of bank position may yield implausible large values
%      if the bottom slope is low or if the GPS coordinates contain an
%      outlier
% TODO not accounting for progressive velocity processing, just using last
%       velocity

    % automatically select all trancsts
    if (nargin() < 2 || isempty(sections))
            sections = 1:length(msh);
    end

    % relative size of cross section that is used to determination the bank position
    %nfrac=[0.1 0.3 0.3 0.5 0.3 0.1 0.1] ;
    if (nargin()<3 || isempty(nfrac))
            nfrac = 0.2*ones(size(sections));
    end

    sfield = 'cs';
   
    % estimate discharge
    for ct=sections(:)'
        [A(ct), p  ] = integrate_area(msh(ct).p, nfrac(ct));
        Q(ct) = integrate_discharge(A(ct), p, msh(ct).(sfield).pars);
        msh(ct).p    = p;
        msh(ct).cs.Q = Q(ct);
        msh(ct).cs.A = A(ct);
    end % for ct
end % function compQ

function [A, p] = integrate_area(p, nfrac)
    % estimate surface area of cells
    p.A=0.5*(sum(bsxfun(@times,p.Z([6 1 5 2],:,end),[1 -1 1 -1]'),1).*...
                    sum(bsxfun(@times,p.N([2 1],:),[1 -1]'),1)+...
                    sum(bsxfun(@times,p.Z([5 2 4 3],:,end),[1 -1 1 -1]'),1).*...
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
    fdx = isnan(p.Atop);
    p.Atop(fdx)     = 0;
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

function Q = integrate_discharge(A,p,vel)
    % Extrapolate velocity and fill in missing velocity
    xvel = vel(p.progfgood_vec(:,end,1));
    fgood  = isfinite(xvel);
    fbad   = ~fgood;
    X      = [p.N(fgood),p.Z(fgood)];
    Nbot = p.Nbot(2,:);
    Zbot = nanmean(p.Zbot([2 5],:),1);
    Ntop = Nbot; % TODO : why not top ?
    Ztop = nanmean(p.Ztop([2 3],:),1);
    Nleft    = nanmean(p.Nleft(1:3,:),1);
    Zleft    = nanmean(p.Zleft(1:3,:),1);
    Nright   = nanmean(p.Nright(1:3,:),1);
    Zright   = nanmean(p.Zright(1:3,:),1);

    % extrapolate top, bottom and side flow
    if exist('scatteredInterpolant','class')==8
        TU = scatteredInterpolant(X,xvel(fgood),'nearest','nearest');
    elseif exist('TriScatteredInterp','class')==8
        TU = TriScatteredInterp(X,xvel(fgood),'nearest'); %#ok<DTRIINT>
    else
        TU = @(x) griddata(X(:,1),X(:,2),xvel(fgood),x(:,1), x(:,2), 'nearest');
    end

    velBot     = TU([Nbot;   Zbot]');
    velTop     = TU([Ntop;   Ztop]');
    velLeft    = TU([Nleft;  Zleft]');
    velRight   = TU([Nright; Zright]');

    % interpolate invalid cells
    xvel(fbad) = TU([p.N(fbad) p.Z(fbad)]);

    % compute discharge
    Q.mid   = sum(xvel.*p.A');
    Q.top   = sum(velTop.*p.Atop');
    Q.bot   = sum(velBot.*p.Abot');
    Q.left  = velLeft.*A.left;
    Q.right = velRight.*A.right;
    Q.total = Q.mid + Q.top + Q.bot + Q.left + Q.right;
end % integrate_discharge

