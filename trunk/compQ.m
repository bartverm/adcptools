clearvars

% load sensitivity
load b_procvel

%% estimate discharge
nfrac=[0.1 0.3 0.3 0.5 0.3 0.1 0.1] ;

for ct=1:7
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
    maxz=nanmax(msh(ct).p.Z(:));
    msh(ct).p.Ztop=repmat([maxz maxz 0 0 maxz]',[1 size(msh(ct).p.Ntop,2)]);
    msh(ct).p.Atop=.5*(sum(bsxfun(@times,msh(ct).p.Ztop([4 1],:),[1 -1]'),1).*...
                       sum(bsxfun(@times,msh(ct).p.Ntop([2 1],:),[1 -1]'),1));
    msh(ct).Ntop=msh(ct).Nbot;
    msh(ct).Ztop=nanmean(msh(ct).p.Ztop([2 3],:),1);
    
    
    % Make cells and determine area for bottom cells
    ntr=round(diff(msh(ct).p.nbed([1 end]))*nfrac(ct)/nanmedian(diff(msh(ct).p.nbed)));
    p=polyfit(msh(ct).p.zbed(1:ntr),msh(ct).p.nbed(1:ntr),1);
    nleft=p(2);
    if ct==5, nleft=-10; end
    p=polyfit(msh(ct).p.zbed(end-ntr+1:end),msh(ct).p.nbed(end-ntr+1:end),1);
    nright=p(2);
    msh(ct).p.Nleft=[nleft; msh(ct).p.nbed([1 1])'; nleft];
    msh(ct).p.Zleft=[0 msh(ct).p.zbed(1) 0 0]';
    msh(ct).Nleft=nanmean(msh(ct).p.Nleft(1:3,:),1);
    msh(ct).Zleft=nanmean(msh(ct).p.Zleft(1:3,:),1);
    msh(ct).p.Aleft=-.5*diff(msh(ct).p.Nleft([1 2]))*msh(ct).p.zbed(1);
    msh(ct).p.Nright=[msh(ct).p.nbed(end) nright msh(ct).p.nbed([end end])]';
    msh(ct).p.Zright=[msh(ct).p.zbed(end) 0 0 msh(ct).p.zbed(end)]';
    msh(ct).p.Aright=.5*diff(msh(ct).p.Nright([1 2]))*diff(msh(ct).p.Zright([1 2]));
    msh(ct).Nright=nanmean(msh(ct).p.Nright(1:3,:),1);
    msh(ct).Zright=nanmean(msh(ct).p.Zright(1:3,:),1);
    
    
    % Extrapolate velocity and fill in missing velocity
    fgood=isfinite(msh(ct).cs.vel(:,:,1));
    fbad=~fgood;
    fgood3=cat(3,fgood,false(size(fgood)),false(size(fgood)));
    fbad3=cat(3,fbad,false(size(fgood)),false(size(fgood)));   
    X=[msh(ct).N(fgood),msh(ct).Z(fgood)];
    TU=scatteredInterpolant(X,msh(ct).cs.vel(fgood3),'nearest','nearest');
    msh(ct).cs.velBot=TU([msh(ct).Nbot;msh(ct).Zbot]');
    msh(ct).cs.velTop=TU([msh(ct).Ntop;msh(ct).Ztop]');
    msh(ct).cs.velLeft=TU([msh(ct).Nleft;msh(ct).Zleft]');
%     fvelleft=any(isfinite(msh(ct).cs.vel
    msh(ct).cs.velRight=TU([msh(ct).Nright;msh(ct).Zright]');
    msh(ct).cs.vel(fbad3)=TU([msh(ct).N(fbad) msh(ct).Z(fbad)]);
    
    % compute Q
    msh(ct).cs.Qmid=sum(msh(ct).cs.vel(msh(ct).p.fgood_3(:,1)).*msh(ct).p.A');
    msh(ct).cs.Qtop=sum(msh(ct).cs.velTop.*msh(ct).p.Atop');
    msh(ct).cs.Qbot=sum(msh(ct).cs.velBot.*msh(ct).p.Abot');
    msh(ct).cs.Qleft=msh(ct).cs.velLeft.*msh(ct).p.Aleft;
    msh(ct).cs.Qright=msh(ct).cs.velRight.*msh(ct).p.Aright;
    msh(ct).cs.Q=msh(ct).cs.Qmid+msh(ct).cs.Qtop+msh(ct).cs.Qbot+msh(ct).cs.Qleft+msh(ct).cs.Qright;
    
    Q(ct)=msh(ct).cs.Q;
    % plot
    figure
    patch(msh(ct).p.N,msh(ct).p.Z,msh(ct).cs.vel(msh(ct).p.fgood_3(:,1))','linestyle','none')
    hold on
    patch(msh(ct).p.Nbot,msh(ct).p.Zbot,msh(ct).cs.velBot')
    patch(msh(ct).p.Ntop,msh(ct).p.Ztop,msh(ct).cs.velTop')
    patch(msh(ct).p.Nleft,msh(ct).p.Zleft,msh(ct).cs.velLeft')
    patch(msh(ct).p.Nright,msh(ct).p.Zright,msh(ct).cs.velRight')
    axis equal
end

