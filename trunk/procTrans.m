function msh=procTrans(adcp,tid,varargin)
% msh=procTrans(ADCP,TID,...)

%    Copyright 2013 Bart Vermeulen
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

% TO DO:
%   - Time variation
%   - Check plotting
%   - Residual removal
%   - Improve depth matching (now a bit slow and possibly inaccurate)
%   - Find direction minizing bathymetry variance?

P=inputParser;
P.FunctionName='procTrans';
P.addParamValue('DepthTransducer',0.3,@(x) isscalar(x) && isnumeric(x) && x>0);
P.addParamValue('DeltaN',5,@(x) isscalar(x) && isnumeric(x) && x>0);
P.addParamValue('DeltaZ',1,@(x) isscalar(x) && isnumeric(x) && x>0);
P.addParamValue('MinimumSigma',0.04,@(x) isscalar(x) && isnumeric(x) && x>=0 && x <=1);
P.parse(varargin{:});

depthtransd=P.Results.DepthTransducer;
veldn=P.Results.DeltaN;
veldz=P.Results.DeltaZ;
sigmin=P.Results.MinimumSigma;

%% Magic numbers
% dprox=20;

%% fetch some adcp data
[x,y]=utmADCP(adcp); % Get adcp position
misal=getExtMisalign(adcp); % Get beam 3 misalignment

% depth location
[dx,dy,dz]=depthADCP(adcp,'Beam3Misalign',misal); % Get offsets from ADCP to measured depths
dz=dz-depthtransd; % transform vertical offset to match with water surface
dx=bsxfun(@plus,x',dx); % transform offsets to depth to UTM coordinate system for horizontal
dy=bsxfun(@plus,y',dy); % transform offsets to depth to UTM coordinate system for horizontal
dxe=nanmean(dx,3);
dye=nanmean(dy,3);
dze=nanmean(dz,3);


% velocity location
[mx,my,mz]=mapADCP(adcp,'IsUpward', false, 'Beam3Misalign',misal); % get offsets to velocity locations
mz=mz-depthtransd; % Correct for depth of transducer
mx=bsxfun(@plus,x',mx); % transform offsets to vel data to UTM coordinate system for horizontal
my=bsxfun(@plus,y',my); % transform offsets to vel data to UTM coordinate system for horizontal
mxe=nanmean(mx,3);
mye=nanmean(my,3);
mze=nanmean(mz,3);


% velocity data
[adcp.VEL,adcp.btvel] = filterADCP2(adcp,'edcf','echotres',20,'difecho',50,'cortres',70,'filterBT',true); % Filter velocity
[vele,btvele]=coradcp(adcp,'e','UseExtHeading',true,'Beam3Misalign',misal); % transform to beam velocity
[adcp.VEL, adcp.btvel]=coradcp(adcp,'b','UseExtHeading',true,'Beam3Misalign',misal); % transform to beam velocity
vel=adcp.VEL-repmat(shiftdim(adcp.btvel,-1),[size(adcp.VEL,1),1,1]); % remove bt from normal vel
vele=vele-repmat(shiftdim(btvele,-1),[size(vele,1),1,1]); % remove bt from normal vel

% tilt information
heading=(getADCPHeading(adcp)+misal)/180*pi; % Get external heading and correct for misalignment
pitch=double(adcp.pitch)/100/180*pi; % Get raw pitch 
roll=double(adcp.roll)/100/180*pi; % Get roll
pitch=atan(tan(pitch).*cos(roll)); % Correct raw pitch to obtain real pitch

bet=20/180*pi;
cb=cos(bet);
sb=sin(bet);
cp=cos(pitch);
sp=sin(pitch);
cr=cos(roll);
sr=sin(roll);
ch=cos(heading');
sh=sin(heading');

TM=cat(3,cat(4,cb.*(ch.*sr - cr.*sh.*sp) + sb.*(ch.*cr + sh.*sp.*sr), - cb.*(sh.*sr + ch.*cr.*sp) - sb.*(cr.*sh - ch.*sp.*sr), cb.*cp.*cr - sb.*cp.*sr),...
         cat(4,cb.*(ch.*sr - cr.*sh.*sp) - sb.*(ch.*cr + sh.*sp.*sr),   sb.*(cr.*sh - ch.*sp.*sr) - cb.*(sh.*sr + ch.*cr.*sp), cb.*cp.*cr + sb.*cp.*sr),...
         cat(4,cb.*(ch.*sr - cr.*sh.*sp) - sb.*cp.*sh, - cb.*(sh.*sr + ch.*cr.*sp) - sb.*ch.*cp, cb.*cp.*cr - sb.*sp),...
         cat(4,cb.*(ch.*sr - cr.*sh.*sp) + sb.*cp.*sh, sb.*ch.*cp - cb.*(sh.*sr + ch.*cr.*sp), sb.*sp + cb.*cp.*cr));

cTM1=repmat(TM(:,:,:,1),[size(vel,1),1,1,1]);
cTM2=repmat(TM(:,:,:,2),[size(vel,1),1,1,1]);
cTM3=repmat(TM(:,:,:,3),[size(vel,1),1,1,1]);

clear heading pitch roll bet cb sb cp sp cr sr ch sh TM

%% find N
msh=repmat(struct(),[size(tid,1),1]); % Initialize final structure for output
n=nan(size(x));
dn=nan(size(dx)); 
mn=nan(size(mx));
ds=nan(size(dx));
ms=nan(size(mx));
dne=nan(size(dxe)); 
mne=nan(size(mxe));
dse=nan(size(dxe));
mse=nan(size(mxe));
T=nan(2,size(tid,1));
N=nan(2,size(tid,1));
Pm=nan(2,size(tid,1));
% Pmid_cell=cell(size(tid,1),1);
for ct=1:size(tid,1)
    fcur=tid(ct,:)>0;
    P=[x(fcur)';y(fcur)'];
    Pm(:,ct)=nanmean(P,2);
    T(:,ct)=princdir(P);
    N(:,ct)=[T(2,ct); -T(1,ct)];
    n(fcur)=T(1,ct)*(x(fcur)-Pm(1,ct))+T(2,ct)*(y(fcur)-Pm(2,ct));
    dn(:,fcur,:)=T(1,ct)*(dx(:,fcur,:)-Pm(1,ct))+T(2,ct)*(dy(:,fcur,:)-Pm(2,ct));
    ds(:,fcur,:)=N(1,ct)*(dx(:,fcur,:)-Pm(1,ct))+N(2,ct)*(dy(:,fcur,:)-Pm(2,ct));
    mn(:,fcur,:)=T(1,ct)*(mx(:,fcur,:)-Pm(1,ct))+T(2,ct)*(my(:,fcur,:)-Pm(2,ct));
    ms(:,fcur,:)=N(1,ct)*(mx(:,fcur,:)-Pm(1,ct))+N(2,ct)*(my(:,fcur,:)-Pm(2,ct));
    dne(:,fcur,:)=T(1,ct)*(dxe(:,fcur,:)-Pm(1,ct))+T(2,ct)*(dye(:,fcur,:)-Pm(2,ct));
    dse(:,fcur,:)=N(1,ct)*(dxe(:,fcur,:)-Pm(1,ct))+N(2,ct)*(dye(:,fcur,:)-Pm(2,ct));
    mne(:,fcur,:)=T(1,ct)*(mxe(:,fcur,:)-Pm(1,ct))+T(2,ct)*(mye(:,fcur,:)-Pm(2,ct));
    mse(:,fcur,:)=N(1,ct)*(mxe(:,fcur,:)-Pm(1,ct))+N(2,ct)*(mye(:,fcur,:)-Pm(2,ct));
    msh(ct).Tvec=T(:,ct);
    msh(ct).Nvec=N(:,ct);
    msh(ct).Pm=Pm(:,ct);
end


%% Get depth averaged velocity
% Calculate sigma for velocity locations
dfg=isfinite(dn) & isfinite(ds) & isfinite(dz);
mfg=isfinite(mn) & isfinite(ms);
md=nan(size(mx));
dfge=isfinite(dne) & isfinite(dse) & isfinite(dze);
mfge=isfinite(mne) & isfinite(mse);
mde=nan(size(mxe));
for ct=1:size(tid,1)
    fcur=tid(ct,:)>0;
    Int=scatteredInterpolant(ds(bsxfun(@and,fcur,dfg)),dn(bsxfun(@and,fcur,dfg)),dz(bsxfun(@and,fcur,dfg)),'natural','linear');
    md(bsxfun(@and,fcur,mfg))=Int(ms(bsxfun(@and,fcur,mfg)),mn(bsxfun(@and,fcur,mfg)));
    Inte=scatteredInterpolant(dse(bsxfun(@and,fcur,dfge))',dne(bsxfun(@and,fcur,dfge))',dze(bsxfun(@and,fcur,dfge))','natural','linear');
    mde(bsxfun(@and,fcur,mfge))=Inte(mse(bsxfun(@and,fcur,mfge)),mne(bsxfun(@and,fcur,mfge)));
end
mSig=1-mz./md;
f_incs=mSig<1 & mSig>sigmin;
% f_close=abs(ms)<dprox;
mSige=1-mze./mde;
f_incse=mSige<1 & mSige>sigmin;
% f_closee=abs(mse)<dprox;

% create n index
ncells=nan(size(tid,1),1);
Pn=nan(2,size(tid,1));
for ct=1:size(tid,1)
    fcur=tid(ct,:)>0;
    nmin=nanmin(mn(bsxfun(@and,fcur,f_incs)));
    Pn(:,ct)=Pm(:,ct)+T(:,ct)*nmin;
    mn(:,fcur,:)=mn(:,fcur,:)-nmin;
    dn(:,fcur,:)=dn(:,fcur,:)-nmin;
    mne(:,fcur,:)=mne(:,fcur,:)-nmin;
    dne(:,fcur,:)=dne(:,fcur,:)-nmin;
    nmax=nanmax(mn(bsxfun(@and,fcur,f_incs)));
    ncells(ct)=ceil(nmax./veldn);
end

nIdx=floor(mn./veldn)+1;
nIdxe=floor(mne./veldn)+1;


% [davel, dastd, damse, daNvel, dax, day]=deal(cell(size(tid,1),1));
% [davele, dastde, daxe, daye]=deal(cell(size(tid,1),1));

% for ct=1:size(tid,1)
%     fcur=bsxfun(@and,f_incs & f_close,tid(ct,:)>0);
%     fcure=bsxfun(@and,f_incse & f_closee,tid(ct,:)>0);
%     dax{ct}=accumarray(nIdx(fcur),mx(fcur),[],@nanmean,nan);
%     day{ct}=accumarray(nIdx(fcur),my(fcur),[],@nanmean,nan);
%     daxe{ct}=accumarray(nIdxe(fcure),mxe(fcure),[],@nanmean,nan);
%     daye{ct}=accumarray(nIdxe(fcure),mye(fcure),[],@nanmean,nan);
%     pIdxe=bsxfun(@plus,nIdxe*0,cat(3,1,2,3));
%     davele{ct}=accumarray({nIdxe(fcure),pIdxe(fcure)},vele(fcure),[],@nanmean,nan);
%     dastde{ct}=accumarray({nIdxe(fcure),pIdxe(fcure)},vele(fcure),[],@nanstd,nan);
%     datacol=cellfun(@(x,y,z,w) [x y z w],...
%         accumarray(nIdx(fcur),vel(fcur),[],@(x) {x},{}),...
%         accumarray(nIdx(fcur),cTM1(fcur),[],@(x) {x},{}),...
%         accumarray(nIdx(fcur),cTM2(fcur),[],@(x) {x},{}),...
%         accumarray(nIdx(fcur),cTM3(fcur),[],@(x) {x},{}),'UniformOutput',false); % Collect beam velocity data and corresponding transformation matrix terms
%     [davel{ct}, dastd{ct}, damse{ct}, daNvel{ct}] =cellfun(@estvel,datacol,'UniformOutput',false); % Estimate velocity with least squares estimates (see function estvel.m for procedure)
%     davel{ct}=squeeze(cell2mat(davel{ct})); % reshape velocity data
% end

%% create mesh
[nbnds, ncntr, d_ncntr, d_nbnds]=deal(cell(size(tid,1),1));
nIdxd=floor((dn+veldn/4)/(veldn/2))+1;
fgood=isfinite(nIdxd) & nIdxd>0 & abs(ds)<10;
fleft=mod(floor(mn/veldn*2),2)==0;
flefte=mod(floor(mne/veldn*2),2)==0;
fright=~fleft;
frighte=~flefte;
dSig=nan(size(mn));
sigIdx=nan(size(mn));
dSige=nan(size(mne));
sigIdxe=nan(size(mne));

for ct=1:size(tid,1)
    fcur=tid(ct,:)>0;
    nbnds{ct}=(0:ncells(ct))*veldn;
    ncntr{ct}=((0:ncells(ct)-1)+0.5)*veldn;
    msh(ct).p.nbed=(0:0.5:ncells(ct))*veldn;
    msh(ct).p.xbed=Pn(1,ct)+T(1,ct)*msh(ct).p.nbed;
    msh(ct).p.ybed=Pn(2,ct)+T(2,ct)*msh(ct).p.nbed;
    dtmp=accumarray(nIdxd(bsxfun(@and,fgood,fcur)),dz(bsxfun(@and,fgood,fcur)),[ncells(ct)*2+1,1],@nanmean,nan);
    fgoodd=isfinite(dtmp);
    Int=griddedInterpolant(find(fgoodd),dtmp(fgoodd));
    dtmp(~fgoodd)=Int(find(~fgoodd)); %#ok<FNDSB>
    msh(ct).p.zbed=dtmp';
    d_nbnds{ct}=dtmp(1:2:end);
    d_ncntr{ct}=dtmp(2:2:end);
    minz_bnd=(1-sigmin)*d_nbnds{ct}; % minimum z for each vertical
    minz_left=minz_bnd(1:end-1);
    minz_right=minz_bnd(2:end);
    minz_cnt=(1-sigmin)*d_ncntr{ct};
    [maxsig, maxloc]=nanmax(mSig(bsxfun(@and,fcur,f_incs)));
    idxf=find(bsxfun(@and,fcur,f_incs));
    maxz=(1-maxsig)*d_ncntr{ct}(nIdx(idxf(maxloc)));
    nz_cnt=round((maxz-minz_cnt)/veldz); % best guess number of z values in vertical giving a dz as close as possible to the given one
    dzed_cnt=(maxz-minz_cnt)./nz_cnt; % compute final dz
    dzed_left=(maxz-minz_left)./nz_cnt;
    dzed_right=(maxz-minz_right)./nz_cnt;
    msh(ct).Z=cumsum([minz_cnt'+dzed_cnt'/2; repmat(dzed_cnt',nanmax(nz_cnt)-1,1)]); % generate z positions
    msh(ct).N=repmat(ncntr{ct},nanmax(nz_cnt),1);
    msh(ct).X=Pn(1,ct)+T(1,ct)*msh(ct).N;
    msh(ct).Y=Pn(2,ct)+T(2,ct)*msh(ct).N;
    msh(ct).Z(bsxfun(@gt,msh(ct).Z,maxz))=nan; % remove points above maxz
   
    % for visualization (compute vertices of mesh cells (previous vectors/matrices give the centers of these cells)
    ZLmin=cumsum([minz_left'; repmat(dzed_left',nanmax(nz_cnt)-1,1)]);
    ZLmax=cumsum([minz_left'+dzed_left'; repmat(dzed_left',nanmax(nz_cnt)-1,1)]);
    ZMmin=cumsum([minz_cnt'; repmat(dzed_cnt',nanmax(nz_cnt)-1,1)]);
    ZMmax=cumsum([minz_cnt'+dzed_cnt'; repmat(dzed_cnt',nanmax(nz_cnt)-1,1)]);
    ZRmin=cumsum([minz_right'; repmat(dzed_right',nanmax(nz_cnt)-1,1)]);
    ZRmax=cumsum([minz_right'+dzed_right'; repmat(dzed_right',nanmax(nz_cnt)-1,1)]);
    LN=repmat(nbnds{ct}(1:end-1),size(ZLmin,1),1); % n of left vertices (same for upper and lower)
    MN=repmat(ncntr{ct},size(ZLmin,1),1); % n of left vertices (same for upper and lower)
    RN=repmat(nbnds{ct}(2:end),size(ZLmin,1),1); % n of right vertices (same for upper and lower)
    nz_cnt(isnan(nz_cnt))=0;
    msh(ct).p.fgood=bsxfun(@le,cumsum(ones(size(msh(ct).Z)),1),nz_cnt');
    msh(ct).p.N=[LN(msh(ct).p.fgood)';...
                 MN(msh(ct).p.fgood)';...
                 RN(msh(ct).p.fgood)';...
                 RN(msh(ct).p.fgood)';...
                 MN(msh(ct).p.fgood)';...
                 LN(msh(ct).p.fgood)';...
                 LN(msh(ct).p.fgood)']; % Generate matrix with n positions as is needed for patch (each column is one cell, each row a vertex)
    msh(ct).p.X=Pn(1,ct)+T(1,ct)*msh(ct).p.N;
    msh(ct).p.Y=Pn(2,ct)+T(2,ct)*msh(ct).p.N;
    msh(ct).p.Z=[ZLmin(msh(ct).p.fgood)';...
                 ZMmin(msh(ct).p.fgood)';...
                 ZRmin(msh(ct).p.fgood)';...
                 ZRmax(msh(ct).p.fgood)';...
                 ZMmax(msh(ct).p.fgood)';...
                 ZLmax(msh(ct).p.fgood)';...
                 ZLmin(msh(ct).p.fgood)']; % Generate matrix with z positions as is needed for patch (each column is one cell, each row a vertex)
     [~, nidx]=ind2sub(size(msh(ct).Z),find(msh(ct).p.fgood));
     msh(ct).p.ZBED=[d_nbnds{ct}(nidx)';...
                     d_ncntr{ct}(nidx)';...
                     d_nbnds{ct}(nidx+1)';...
                     d_nbnds{ct}(nidx+1)';...
                     d_ncntr{ct}(nidx)';...
                     d_nbnds{ct}(nidx)';...
                     d_nbnds{ct}(nidx)'];                   
     msh(ct).p.Sig=1-msh(ct).p.Z./msh(ct).p.ZBED;
%     nel2=numel(msh(ct).p.zbed)/2; % compute element idx halfway the nvector
%     msh(ct).p.n_watsurf=[interp1(msh(ct).p.zbed(1:floor(nel2)),msh(ct).p.nvec(1:floor(nel2)),0),...
%                      interp1(msh(ct).p.zbed(ceil(nel2):end),msh(ct).p.nvec(ceil(nel2):end),0)]; % Interpolate left and right n-coordinate of water surface line
    dSigl=-dzed_left./d_nbnds{ct}(1:end-1);
    dSigc=-dzed_cnt./d_ncntr{ct}(1:end);
    dSigr=-dzed_right./d_nbnds{ct}(2:end);

    fnd=bsxfun(@and,tid(ct,:)>0,f_incs & fleft);
    dSig(fnd)=dSigl(nIdx(fnd))+(dSigc(nIdx(fnd))-dSigl(nIdx(fnd)))/veldn*2.*(mn(fnd)-nbnds{ct}(nIdx(fnd))');
    fnd=bsxfun(@and,tid(ct,:)>0,f_incs & fright);
    dSig(fnd)=dSigc(nIdx(fnd))+(dSigr(nIdx(fnd))-dSigc(nIdx(fnd)))/veldn*2.*(mn(fnd)-ncntr{ct}(nIdx(fnd))');

    fnde=bsxfun(@and,tid(ct,:)>0,f_incse & flefte);
    dSige(fnde)=dSigl(nIdxe(fnde))+(dSigc(nIdxe(fnde))-dSigl(nIdxe(fnde)))/veldn*2.*(mne(fnde)-nbnds{ct}(nIdxe(fnde))');
    fnde=bsxfun(@and,tid(ct,:)>0,f_incse & frighte);
    dSige(fnde)=dSigc(nIdxe(fnde))+(dSigr(nIdxe(fnde))-dSigc(nIdxe(fnde)))/veldn*2.*(mne(fnde)-ncntr{ct}(nIdxe(fnde))');

    fnd=bsxfun(@and,tid(ct,:)>0,f_incs);
    sigIdx(fnd)=floor((mSig(fnd)-sigmin)./dSig(fnd))+1;
    
    fnde=bsxfun(@and,tid(ct,:)>0,f_incse);
    sigIdxe(fnde)=floor((mSige(fnde)-sigmin)./dSige(fnde))+1;

    % velocity
    
    fnd=bsxfun(@and,tid(ct,:)>0,f_incs & sigIdx<=size(msh(ct).Z,1));
%     msh(ct).Xav=accumarray({sigIdx(fg) betIdx(fg)},cx(fg),size(msh(ct).Sig),@nanmean,nan);
%     msh(ct).Yav=accumarray({sigIdx(fg) betIdx(fg)},cy(fg),size(msh(ct).Sig),@nanmean,nan);
%     msh(ct).Zav=accumarray({sigIdx(fg) betIdx(fg)},cz(fg),size(msh(ct).Sig),@nanmean,nan);
    msh(ct).datacol=cellfun(@(x,y,z,w) [x y z w],...
        accumarray({sigIdx(fnd) nIdx(fnd)},vel(fnd),size(msh(ct).Z),@(x) {x},{}),...
        accumarray({sigIdx(fnd) nIdx(fnd)},cTM1(fnd),size(msh(ct).Z),@(x) {x},{}),...
        accumarray({sigIdx(fnd) nIdx(fnd)},cTM2(fnd),size(msh(ct).Z),@(x) {x},{}),...
        accumarray({sigIdx(fnd) nIdx(fnd)},cTM3(fnd),size(msh(ct).Z),@(x) {x},{}),'UniformOutput',false); % Collect beam velocity data and corresponding transformation matrix terms
    
    [msh(ct).vel, msh(ct).std, msh(ct).mse, msh(ct).Nvel] =cellfun(@estvel,msh(ct).datacol,'UniformOutput',false); % Estimate velocity with least squares estimates (see function estvel.m for procedure)
%     
    msh(ct).vel=squeeze(cell2mat(msh(ct).vel)); % reshape velocity data
    msh(ct).mse=squeeze(cell2mat(msh(ct).mse)); % reshape mean squared error data
    msh(ct).std=squeeze(cell2mat(msh(ct).std)); % reshape standard deviation data
    msh(ct).Nvel=squeeze(cell2mat(msh(ct).Nvel)); % reshape Number of valid samples data

    fnde=bsxfun(@and,tid(ct,:)>0,repmat(f_incse & sigIdxe<=size(msh(ct).Z,1),[1 1 3]));
    tsigIdxe=repmat(sigIdxe,[1 1 3]);
    tnIdxe=repmat(nIdxe,[1 1 3]);
    pIdxe=bsxfun(@plus,tnIdxe*0,cat(3,1,2,3));
    msh(ct).vele=accumarray({tsigIdxe(fnde), tnIdxe(fnde), pIdxe(fnde)},vele(fnde),[size(msh(ct).Z) 3],@nanmean,nan);
    msh(ct).stde=accumarray({tsigIdxe(fnde), tnIdxe(fnde), pIdxe(fnde)},vele(fnde),[size(msh(ct).Z) 3],@nanstd,nan);
    msh(ct).Nvele=accumarray({tsigIdxe(fnde), tnIdxe(fnde), pIdxe(fnde)},vele(fnde),[size(msh(ct).Z) 3],@numel,0);
end

function [vel, std, mse, nin]=estvel(in)
% estimates cartesian velocities given beam velocities and correspoding
% transformation matrix terms. IN is Nx4, N being the amount of beam
% velocity samples. Each row is composed of a beam velocity and the three
% terms from the transformation matrix

in(any(isnan(in),2),:)=[]; % Throw out bad data
%rmvel=true; % initialize variable holding indices to outliers

% while any(rmvel) % iterate until no outliers are found
    if isempty(in) || rank(in(:,2:4),1)<3 % if no velocity is available or matrix is rank deficient
        vel=[nan;nan;nan];std=vel; mse=nan; % return nans
%       break; % We're done here
    else
        [vel, std, mse]=lscov(in(:,2:4),in(:,1)); % Perform least squares estimate of velocities
    end
%    res=abs(in(:,1)-in(:,2:4)*vel); % Determine residuals in beam velocity
%    res>6*median(res); % Detect outliers from residuals
%     in(rmvel,:)=[]; % Remove outliers
%end
nin=shiftdim(size(in,1),-3); % compute amount of points used for final velocity computation
vel=shiftdim(vel,-3); % shift dimension to output 1x1x1x3 vector
std=shiftdim(std,-3); % shift dimension to output 1x1x1x3 vector
mse=shiftdim(mse,-3); % shift dimension to output 1x1x1x1 scalar



function [V, mag]=princdir(P,varargin)
% PRINCDIR find the principal vector of a scattered cloud of points. 
%   V=PRINCDIR(P) returns the vector V pointing in the direction of the 
%       largest variance for the cloud of points defined in P. P is a 
%       matrix with one point per column and the amount of rows 
%       corresponding to the amount of dimensions.
%
%   [V, var]=PRINCDIR(P) also return the variance in the direction of V
%
%   PRINCDIR(P,dim) indicates the dimension DIM along which the vectors are
%       defined.
%
%   Author: Bart Vermeulen
%   Last edit: 07-08-2013


assert(isnumeric(P) && ismatrix(P));

dim=1;
if nargin==2
    dim=varargin{1};
    assert(isequal(dim,1) || isequal(dim,2))
end

fg=all(isfinite(P),dim);

if dim==1
    C=cov(P(:,fg)');
else
    C=cov(P(fg,:));
end

[eigvec, eigval]=eig(C);

[mag, fprinc]=max(sum(eigval,2));

V=eigvec(:,fprinc);
