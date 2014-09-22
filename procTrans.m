function msh=procTrans(adcp,tid,varargin)
% procTrans processes ADCP transect data
%   MSH=procTrans(ADCP,TID,...) process the data in the ADCP structure (as
%       read by readADCP. TID identifies sections and crossings. It is a
%       NxM matrix with N equal to the amount of sections and M equal to
%       the amount of ensembles in the ADCP structure. Every element in TID
%       higher than zero at location (n,m) indicates that ensemble m
%       belongs to section n. Each element in TID can contain a number
%       which indicates the crossing an element belongs to.
%       The output MSH is an N by 1 structure containing the results for
%       each section. It contains the following fields:
%           Tvec: The unit vector indicating the tangential direction of 
%               the section
%           Nvec: The unit vector orthogonal to Tvec
%           Pm: Position vector indicating the center of the section
%           eta: Vector containing the water level eta for each time
%               processing step
%           time: Average time of each time processing step
%           Sig,N,X,Y: Coordinates of the centres of the processing cells.
%               These matrices all have dimensions [Nz,Nn], being the
%               maximum number of cells in depth and the number of cells in
%               the width direction. They contain the Sigma coordinates,
%               i.e. the relative elevation, N coordinate parallel to Tvec,
%               X coordinate in UTM and Y coordinate in UTM.
%           Z,T: Vertical coordinate and time for the cell centres. These
%               matrices have dimensions [Nz,Nn,Nt], which includes and
%               additional dimension for the time variation.
%           vel: The velocity estimates. The matrix has dimensions
%               [Nz,Nn,Nt,Nd], i.e. the cells in depth, width, time, and
%               the Nd=3, being the three components of the velocity
%               vector: eastward, northward, upward.
%           std: Standard deviation in the velocity estimates. Has the same
%               dimensions as vel
%           mse: Mean square error of the velocity least square fitting
%           Nvel: Number of individual beam velocity samples used in the
%               velocity estimate
%           vele,stde,Nvele: The same as vel, std and Nvel, but now
%               computed directly combining beam velocity samples, instead
%               of doing this afterwards.
%           prog[vel, std, mse, Nvel, vele, stde, Nvele]: Same as above but
%               now with the processing result progressively including a
%               crossing. The have an additional dimension Nc, equal to the
%               number of crossings. These are only generated when the
%               'Progressive' paramter is set to true.
%       Additionaly the structure contain another structure 'p' which is
%       mainly meant for plotting. This structure contains:
%           nbed, xbed, ybed, zbed: n,x,y and z coordinate of the bed
%               computed at the cell boundaries and centers.
%           fgood: mesh cells that are actually part of the section (not
%               below the bottom or above the surface).
%           N,X,Y,Sig: Coordinates of patches defining mesh cells in the
%               cross-section (indicated by fgood). Size of the matrices is
%               7xNcells, with Ncells being the amount on cells in the
%               section.
%           Z: Vertical coordinate of the patches defining the mesh-cells.
%               This matrix has dimensions 7xNcellsxNt, with Nt being the
%               number of time steps used in the processing
%           
%    MSH=procTrans(ADCP, TID,'PropertyName',PropertyValue) 
%         allows to specify the following additional settings:
%
%         DepthTransducer
%         A numerical scalar indicating the depth of the transducer (in m).
%         Default is 0.3 m.
%
%         DeltaN
%         Cross-section resolution of the mesh (in m). Default is 5 m.
%
%         DeltaZ
%         Vertical resolution of the mesh (in m). Default is 1 m.
%
%         DeltaT
%         Spatial resolution (in minutes). Setting this to 0 indicates
%         temporal averaging over all data (No time variation in the data).
%         Default is 0 minutes.
%
%         Eta
%         Vector with as many elements as ensembles indicating the water
%         level fluctuation in m. Default is 0 (no water level
%         fluctuation).
%
%         MinimumSigma
%         Minimum relative depth to include the data in the processing.
%         Numerical scalar with values ranging from 0 to 1 (although
%         usually close to 0). Default is 0.04.
%
%         EstimateGradients
%         true|{false}
%         Adds spatial and temporal gradients to the fitted model
%
%         Progressive
%         true|{false}
%         Indicates whether the velocity should be solved progressively
%         including more crossings. This is especially usefull to evaluate
%         the quality turbulence averaging.
%
%         RemoveOutliers
%         Numerical scalar Nr, indicating to remove beam velocities which 
%         have a residual exceeding Nr times the median residual of a 
%         velocity estimate. If Nr=0 no outliers are removed. Default is 0.
%
%         StdFiltering
%         Indicates the how many times the standard deviation of a velocity
%         should exceed the median standard deviation to be regarded as
%         bad. Default is 6.
%
%         Proximity
%         Numerical scalar indicating a maximum distance from the section
%         for velocity data to be included in the velocity estimates
%
%         ShipReference
%         {'bt'}|'gps'|'btgps'|'gpsbt'
%         Select reference to use to determine ship velocity: 'bt' for
%         bottom tracking, 'gps' for gps, 'btgps' for bottom tracking,
%         using gps when bottom tracking is missing and 'gpsbt' to use gps
%         and bottom tracking when gps is missing.
%
%   Example: msh=procTrans(ADCP,tid,'Eta',eta,'Progressive',true,'DeltaT',
%         60)
%         Processes the data in ADCP for sections defined in tid. It will
%         process the data progressively including crossing, defined by the
%         numbers in tid. Water level fluctuation is given by eta. Time
%         step will be 60 minutes.
%   
%   Author: Bart Vermeulen


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
%   - Check plotting
%   - Residual removal
%   - Improve depth matching (now a bit slow and possibly inaccurate)
%   - Find direction minizing bathymetry variance?
%   - Detect Eta from bathymetry
%   - Add code documentation
%   - Change to covariance matrix instead of standard deviation
%   - Somewhere p.zbed is generated above water surface... could be a
%       problem!

%% Handle input
P=inputParser;
P.FunctionName='procTrans';
P.addParamValue('DepthTransducer',0.3,@(x) isscalar(x) && isnumeric(x) && x>0);
P.addParamValue('DeltaN',5,@(x) isscalar(x) && isnumeric(x) && x>0);
P.addParamValue('DeltaZ',1,@(x) isscalar(x) && isnumeric(x) && x>0);
P.addParamValue('DeltaT',0,@(x) isscalar(x) && isnumeric(x) && x>0);
P.addParamValue('Eta',0,@(x) isvector(x) && isnumeric(x))
P.addParamValue('MinimumSigma',0.04,@(x) isscalar(x) && isnumeric(x) && x>=0 && x <=1);
P.addParamValue('Progressive',false,@(x) isscalar(x) && islogical(x));
P.addParamValue('RemoveOutliers',0,@(x) isscalar(x) && isnumeric(x) && x>=0);
P.addParamValue('Proximity',0,@(x) isscalar(x) && isnumeric(x) && x>=0);
P.addParamValue('StdFiltering',6,@(x) isscalar(x) && isnumeric(x) && x>=0);
P.addParamValue('ShipReference','bt',@(x) ischar(x) && any(strcmpi(x,{'bt','gps','gpsbt','btgps'})));
P.addParamValue('Pusr',[], @(x) isnumeric(x));
P.parse(varargin{:});

assert(isADCPstruct(adcp));
assert(size(tid,2)==size(adcp.VEL,2));

depthtransd=P.Results.DepthTransducer;
veldn=P.Results.DeltaN;
veldz=P.Results.DeltaZ;
veldt=P.Results.DeltaT/24/60;
eta=P.Results.Eta;
sigmin=P.Results.MinimumSigma;
progflag=P.Results.Progressive;
rmresid=P.Results.RemoveOutliers;
proxim=P.Results.Proximity;
nstd=P.Results.StdFiltering;
time=datenum(adcp.timeV)'; % Get time
eta=eta(:)';    % Vectorize eta
eta=eta.*ones(1,size(numel(time),1)); % Ensure eta has same number of elements as the number of ensembles
shref=P.Results.ShipReference;
Pusr=P.Results.Pusr;

%% Get ADCP positioning data
[x,y]=utmADCP(adcp); % Get adcp position
misal=getExtMisalign(adcp); % Get beam 3 misalignment

% depth location
[dx,dy,dz]=depthADCP(adcp,'Beam3Misalign',misal); % Get offsets from ADCP to measured depths
dz=bsxfun(@plus,dz-depthtransd,eta); % transform vertical offset to match with water surface
dx=bsxfun(@plus,x',dx); % transform offsets to depth to UTM coordinate system for horizontal
dy=bsxfun(@plus,y',dy); % transform offsets to depth to UTM coordinate system for horizontal
dxe=nanmean(dx,3); % determine beam average for conventional velocity processing
dye=nanmean(dy,3); % determine beam average for conventional velocity processing
dze=nanmean(dz,3); % determine beam average for conventional velocity processing


% velocity location
[mx,my,mz]=mapADCP(adcp,'IsUpward', false, 'Beam3Misalign',misal); % get offsets to velocity locations
mz=bsxfun(@plus,mz-depthtransd,eta); % Correct for depth of transducer and eta
mx=bsxfun(@plus,x',mx); % transform offsets to vel data to UTM coordinate system for horizontal
my=bsxfun(@plus,y',my); % transform offsets to vel data to UTM coordinate system for horizontal
mxe=nanmean(mx,3); % determine beam average for conventional velocity processing
mye=nanmean(my,3); % determine beam average for conventional velocity processing
mze=nanmean(mz,3); % determine beam average for conventional velocity processing


% velocity data
[adcp.VEL,adcp.btvel] = filterADCP(adcp,'','filterBT',true); % Filter velocity
gpsvel=getGPSvel(adcp,'b'); 
gpsvele=getGPSvel(adcp,'e');

[vele,btvele]=corADCP(adcp,'e','UseExtHeading',true,'Beam3Misalign',misal); % transform to beam velocity
[adcp.VEL, adcp.btvel]=corADCP(adcp,'b','UseExtHeading',true,'Beam3Misalign',misal); % transform to beam velocity

fbad_bt=any(isnan(adcp.btvel),2);
fbad_bte=any(isnan(btvele),2);
dz(:,fbad_bt,:)=nan;
dze(:,fbad_bte,:)=nan;

if strcmpi(shref,'gps')
    adcp.btvel=gpsvel;
    adcp.btvele=gpsvele;
elseif strcmpi(shref,'btgps')
    adcp.btvel(fbad_bt,:)=gpsvel(fbad_bt,:);
    btvele(fbad_bte,:)=gpsvele(fbad_bte,:);
elseif strcmpi(shref,'gpsbt')
    fbad_gps=any(isnan(gpsvel),2);
    gpsvel(fbad_gps,:)=adcp.btvel(fbad_gps,:);
    adcp.btvel=gpsvel;
    fbad_gpse=any(isnan(gpsvele),2);
    gpsvele(fbad_gpse,:)=btvele(fbad_gpse,:);
    btvele=gpsvele;
end

vel=adcp.VEL-repmat(shiftdim(adcp.btvel,-1),[size(adcp.VEL,1),1,1]); % remove bt from normal vel
vele=vele-repmat(shiftdim(btvele,-1),[size(vele,1),1,1]); % remove bt from normal vel

% tilt information
heading=(getADCPHeading(adcp)+misal)/180*pi; % Get external heading and correct for misalignment
pitch=double(adcp.pitch)/100/180*pi; % Get raw pitch 
roll=double(adcp.roll)/100/180*pi; % Get roll
pitch=atan(tan(pitch).*cos(roll)); % Correct raw pitch to obtain real pitch

bet=20/180*pi; % beam angle = beta
cb=cos(bet); % cos beta
sb=sin(bet); % sin beta
cp=cos(pitch); % cos pitch
sp=sin(pitch); % sin pitch
cr=cos(roll); % cos roll
sr=sin(roll); % sin roll
ch=cos(heading); % cos heading
sh=sin(heading); % sin heading

% Determine transformation matrix (from beam coords to earth coords)
TM=cat(3,cat(4,cb.*(ch.*sr - cr.*sh.*sp) + sb.*(ch.*cr + sh.*sp.*sr), - cb.*(sh.*sr + ch.*cr.*sp) - sb.*(cr.*sh - ch.*sp.*sr), cb.*cp.*cr - sb.*cp.*sr),...
         cat(4,cb.*(ch.*sr - cr.*sh.*sp) - sb.*(ch.*cr + sh.*sp.*sr),   sb.*(cr.*sh - ch.*sp.*sr) - cb.*(sh.*sr + ch.*cr.*sp), cb.*cp.*cr + sb.*cp.*sr),...
         cat(4,cb.*(ch.*sr - cr.*sh.*sp) - sb.*cp.*sh, - cb.*(sh.*sr + ch.*cr.*sp) - sb.*ch.*cp, cb.*cp.*cr - sb.*sp),...
         cat(4,cb.*(ch.*sr - cr.*sh.*sp) + sb.*cp.*sh, sb.*ch.*cp - cb.*(sh.*sr + ch.*cr.*sp), sb.*sp + cb.*cp.*cr));

% Split in three matrices each with one term of the transformation matrix
TM1=repmat(TM(:,:,:,1),[size(vel,1),1,1,1]);
TM2=repmat(TM(:,:,:,2),[size(vel,1),1,1,1]);
TM3=repmat(TM(:,:,:,3),[size(vel,1),1,1,1]);

%cleanup
clear heading pitch roll bet cb sb cp sp cr sr ch sh TM



%% Project positions to section
% Initialize
nsec=size(tid,1);
msh=repmat(struct(),[nsec,1]); % Initialize final structure for output
T=nan(2,size(tid,1)); % Unit vector tangential to section
N=nan(2,size(tid,1)); % Unit vector orthogonal to section 
Pm=nan(2,size(tid,1)); % Average section location vector
Pn=nan(2,size(tid,1)); % Average section location vector


%% create time index (mapping to right time cell)
% Initialize variables
tIdx=ones(size(time));

% Compute time index
if veldt>0 % Only compute if positive time step is given
    tIdx=floor((time-min(time))./veldt)+1;
end

% Compute average time and eta for time cells
time_cnt=shiftdim(accumarray(tIdx',time',[],@nanmean,nan),-2);
eta_cnt=shiftdim(accumarray(tIdx',eta',[],@nanmean,nan),-2);

% Replicate time index matrix
tIdx=repmat(tIdx,[size(mx,1),1,size(mx,3)]);


% Start processing
for ct=1:size(tid,1) % For all sections  
    msh(ct).eta=eta_cnt; % store average eta for time steps
    msh(ct).time=time_cnt; % store average time for time steps

    fcur=tid(ct,:)>0; % mask, indicating data belonging to current section
    
    % Compute relevant vectors
    P=[x(fcur)';y(fcur)']; % ADCP position vectors
    if (~isempty(Pusr))
	% todo, multiple transects
	Pm = mean(Pusr(:,:,ct),2);
	T(:,ct) = diff(Pusr(:,:,ct),[],2);
	T(:,ct) = T(:,ct)/norm(T(:,ct));
    else
    Pm(:,ct)=nanmean(P,2); % Average ADCP position vector
    T(:,ct)=princdir(P); % Section tangential vector, determined as largest eigenvector of P
    end
    N(:,ct)=[T(2,ct); -T(1,ct)]; % Section orthogonal vector (orthogonal to T)
    
    % Project positions (inner product with unit vectors)
%     n=T(1,ct)*(x(fcur)-Pm(1,ct))+T(2,ct)*(y(fcur)-Pm(2,ct)); % ADCP position vector, n component
    dn=T(1,ct)*(dx(:,fcur,:)-Pm(1,ct))+T(2,ct)*(dy(:,fcur,:)-Pm(2,ct)); % depth position, n component
    ds=N(1,ct)*(dx(:,fcur,:)-Pm(1,ct))+N(2,ct)*(dy(:,fcur,:)-Pm(2,ct)); % depth position, s component
    mn=T(1,ct)*(mx(:,fcur,:)-Pm(1,ct))+T(2,ct)*(my(:,fcur,:)-Pm(2,ct)); % velocity position, n component
    ms=N(1,ct)*(mx(:,fcur,:)-Pm(1,ct))+N(2,ct)*(my(:,fcur,:)-Pm(2,ct)); % velocity position, s component
    dne=T(1,ct)*(dxe(:,fcur,:)-Pm(1,ct))+T(2,ct)*(dye(:,fcur,:)-Pm(2,ct)); % depth position, n component (conventional processing)
    dse=N(1,ct)*(dxe(:,fcur,:)-Pm(1,ct))+N(2,ct)*(dye(:,fcur,:)-Pm(2,ct)); % depth position, s component (conventional processing)
    mne=T(1,ct)*(mxe(:,fcur,:)-Pm(1,ct))+T(2,ct)*(mye(:,fcur,:)-Pm(2,ct)); % velocity position, n component (conventional processing)
    mse=N(1,ct)*(mxe(:,fcur,:)-Pm(1,ct))+N(2,ct)*(mye(:,fcur,:)-Pm(2,ct)); % velocity position, s component (conventional processing)
    
    % Store important vectors
    msh(ct).Tvec=T(:,ct); % tangential vector
    msh(ct).Nvec=N(:,ct); % orthogonal vector
    msh(ct).Pm=Pm(:,ct); % mean position vector
    
    %% Calculate sigma for velocity locations
    cdz=dz(:,fcur,:);
    cdze=dze(:,fcur);
    
    % Finite depth and velocity masks
    dfg=isfinite(dn) & isfinite(ds) & isfinite(cdz); % Mask for finite depth locations
    mfg=isfinite(mn) & isfinite(ms); % Mask for finite velocity locations
    dfge=isfinite(dne) & isfinite(dse) & isfinite(cdze); % Mask for finite depth locations (conventional processing)
    mfge=isfinite(mne) & isfinite(mse); % Mask for finite velocity locations (conventional processing)

    % Initialize variables
    md=nan(size(mn)); % Bed elevation at velocity locations
    mde=nan(size(mne)); % Bed elevation at velocity locations (conventional processing)

    % Compute Bed elevations at velocity sample locations (use accumarray instead?)
    if verLessThan('matlab','7.10')
        md(mfg) = griddata(ds(dfg),dn(dfg),cdz(dfg),ms(mfg),mn(mfg),'natural'); %#ok<GRIDD>
        mde(mfg) = griddata(dse(dfge),dne(dfge),cdze(dfge),mse(mfge),mn(mfge),'natural'); %#ok<GRIDD>
    elseif verLessThan('matlab','8.1')
        Int=TriScatteredInterp(ds(dfg),dn(dfg),cdz(dfg),'natural'); %#ok<REMFF1> % Create interpolant (natural interpolation with linear extrapolation)
        md(mfg)=Int(ms(mfg),mn(mfg)); % interpolate bed elevation at velocity locations
        Inte=TriScatteredInterp(dse(dfge)',dne(dfge)',cdze(dfge)','natural'); %#ok<REMFF1> % Create interpolant (natural interpolation with linear extrapolation) (conventional processing)
        mde(mfge)=Inte(mse(mfge),mne(mfge)); % interpolate bed elevation at velocity locations (conventional processing)
   else
        Int=scatteredInterpolant(ds(dfg),dn(dfg),cdz(dfg),'natural','linear'); % Create interpolant (natural interpolation with linear extrapolation)
        md(mfg)=Int(ms(mfg),mn(mfg)); % interpolate bed elevation at velocity locations
        Inte=scatteredInterpolant(dse(dfge)',dne(dfge)',cdze(dfge)','natural','linear'); % Create interpolant (natural interpolation with linear extrapolation) (conventional processing)
        mde(mfge)=Inte(mse(mfge),mne(mfge)); % interpolate bed elevation at velocity locations (conventional processing)
    end  

    % Compute sigma
    cmz=mz(:,fcur,:);
    cmze=mze(:,fcur,:);
    if numel(eta)>1, ceta=eta(fcur); else ceta=eta; end
    mSig=1-bsxfun(@minus,cmz,ceta)./bsxfun(@minus,md,ceta); % Compute sigma for velocity locations
    f_incs=mSig<1 & mSig>sigmin; % Create mask for velocity locations in the cross-section
    mSige=1-bsxfun(@minus,cmze,eta)./bsxfun(@minus,mde,eta); % Compute sigma for velocity locations (conventional processing)
    f_incse=mSige<1 & mSige>sigmin; % Create mask for velocity locations in the cross-section (conventional processing)
    
    %% create n index (maps velocity to correct n vertical)
    if (~isempty(Pusr))
        nmin = min(T(:,ct)'*(Pusr(:,:,ct) - Pm(:,ct)*[1 1]));
    else
    nmin=nanmin(mn(f_incs)); % determine minimum n
    end
    Pn(:,ct)=Pm(:,ct)+T(:,ct)*nmin; % compute vector pointing to minimum n location (projected on section)
    mn=mn-nmin; % remove minimum n from velocity n coordinates (now they range between 0 and maximum n)
    dn=dn-nmin; % remove minimum n from depth n coordinates (now they range between 0 and maximum n)
    mne=mne-nmin; % remove minimum n from velocity n coordinates (now they range between 0 and maximum n) (conventional processing)
%     dne=dne-nmin; % remove minimum n from depth n coordinates (now they range between 0 and maximum n) (conventional processing)
    if (~isempty(Pusr))
        nmax = max(T(:,ct)'*(Pusr(:,:,ct) - Pn(:,ct)*[1 1]));
    else
    nmax=nanmax(mn(f_incs)); % compute maximum n
    end
    ncells=ceil(nmax./veldn); % determine number of vertical in section

    % Compute n indices
    nIdx=floor(mn./veldn)+1; % Compute N index, mapping velocity data to correct vertical
    nIdxe=floor(mne./veldn)+1; % Compute N index, mapping velocity data to correct vertical (conventional processing)

    % select points within transect range
    inrange = nIdx >0 & nIdx <= ncells;
    inrange_e = nIdxe >0 & nIdxe <= ncells;
    fdx = find(nIdx < 1 & nIdx > ncells);
    nIdx(fdx) = NaN;
    fdx = find(nIdxe < 1 & nIdxe > ncells);
    nIdxe(fdx) = NaN;

    %% Meshing
    % Initialize variables
%     [nbnds, ncntr, d_ncntr, d_nleft, d_nright, d_nbnds]=deal(cell(size(tid,1),1)); % Initialize, n coordinate for cell boundaries and center, depth at cell center, at cell boundaries, at cell left boundary and cell right boundary
    dSig=nan(size(mn)); % Step in sigma
    sigIdx=nan(size(mn)); % Sigma index (mapping to right depth cell)
    dSige=nan(size(mne)); % Step in sigma (conventional processing)
    sigIdxe=nan(size(mne)); % Sigma index (mapping to right depth cell) (conventional processing)

    % Preliminary computations
    nIdxd=floor((dn+veldn/4)/(veldn/2))+1; % Index to determine depth at cell boundaries and centres
    fgood=isfinite(nIdxd) & nIdxd>0;% & abs(ds)<2*veldn; % MAGIC NUMBER!!!! % Select data to include in depth calculation at cell edges and center
    fleft=mod(floor(mn/veldn*2),2)==0; % Selects data on the left half of a cell
    flefte=mod(floor(mne/veldn*2),2)==0; % Selects data on the left half of a cell (conventional processing)
    fright=~fleft; % Selects data on the left half of a cell
    frighte=~flefte; % Selects data on the left half of a cell (conventional processing)

    % determine n coordinate and depth of cell edges and centres
    nbnds=(0:ncells)*veldn; % n coordinate of cell boundaries
    ncntr=((0:ncells-1)+0.5)*veldn; % n coordinate of cell centres
    dtmp=accumarray(nIdxd(fgood),cdz(fgood),[ncells*2+1,1],@nanmean,nan); % Determine elevation (z) at cell centres and boundaries (where available)
    fgoodd=isfinite(dtmp); % find all available elevations at cell centres and boundaries
    Int=griddedInterpolant(find(fgoodd),dtmp(fgoodd)); % create interpolant to calculate elevation where it's missing (probably at the edges)
    dtmp(~fgoodd)=Int(find(~fgoodd)); %#ok<FNDSB> % interpolate at cell centres or boundaries with missing depth
    d_ncntr=dtmp(2:2:end); % get elevation at cell centres

    % Compute minimum and maximum z and sigma for verticals
    [maxsig, maxloc]=nanmax(mSig(f_incs)); % maximum Sigma measured in all data belonging to current section
    idxf=find(f_incs); % get indices to which index of maximum computed in line above refers to
    maxz=(1-maxsig)*d_ncntr(nIdx(idxf(maxloc))); % Transform maximum sigma to maximum mesh elevation at section location
    
    % check for degenerate cells, ie. left or right boundary exceeds
    % maximum z
    ntmp=(0:0.5:ncells)*veldn; 
    if (1-sigmin)*dtmp(1)>maxz
        ntmp(1)=ntmp(2)-abs(((1-sigmin)*dtmp(2)-maxz)./(dtmp(2)-dtmp(1))*(veldn/2));
        dtmp(1)=maxz/(1-sigmin);
    end
    if (1-sigmin)*dtmp(end)>maxz
        ntmp(end)=ntmp(end-1)+abs(((1-sigmin)*dtmp(end-1)-maxz)./(dtmp(end)-dtmp(end-1))*(veldn/2));
        dtmp(end)=maxz/(1-sigmin);
    end
    
    msh(ct).p.nbed=ntmp; % n coordinate of both cell centres and boundaries
    msh(ct).p.xbed=Pn(1,ct)+T(1,ct)*msh(ct).p.nbed; % x coordinate of both cell centres and boundaries
    msh(ct).p.ybed=Pn(2,ct)+T(2,ct)*msh(ct).p.nbed; % y coordinate of both cell centres and boundaries
    msh(ct).p.zbed=dtmp'; % store the z coordinate of the bed for cell centres and boundaries

    d_nbnds=dtmp(1:2:end); % get average bed elevation at cell boundaries
    d_nleft=d_nbnds(1:end-1); % get average bed elevation at cell left boundary
    d_nright=d_nbnds(2:end); % get average bed elevation at cell left boundary
    minz_bnd=(1-sigmin)*d_nbnds; % minimum z for each vertical at cell boundaries
    minz_left=minz_bnd(1:end-1); % minimum z for each vertical at cell left boundaries
    minz_right=minz_bnd(2:end); % minimum z for each vertical at cell right boundaries
    minz_cnt=(1-sigmin)*d_ncntr; % minimum z for each vertical at cell centres

    
    maxsig=1-maxz./d_ncntr';

    % Compute number of cells in each vertical and cell-size (in sigma and zed)
    nz_cnt=round((maxz-minz_cnt)/veldz); % best guess number of z values in vertical giving a dz as close as possible to the given one
    dzed_cnt=(maxz-minz_cnt)./nz_cnt; % compute dz at cell centers
    dsig_cnt=-dzed_cnt./d_ncntr; % compute dsigma at cell centers
    dzed_left=(maxz-minz_left)./nz_cnt; % compute dz at left cell boundary
    dsig_left=-dzed_left./d_nleft; % compute dsigma at left cell boundary
    dzed_right=(maxz-minz_right)./nz_cnt; % compute dz at right cell boundary
    dsig_right=-dzed_right./d_nright; % comput dsigma at right cell boundary

    % Generate mesh (only cell centres)
    msh(ct).Sig=cumsum([sigmin+dsig_cnt'/2; repmat(dsig_cnt',nanmax(nz_cnt)-1,1)]); % generate sigma positions of cell centres
    msh(ct).Sig(bsxfun(@gt,msh(ct).Sig,maxsig))=nan; % Remove cell centers outside the section
    msh(ct).N=repmat(ncntr,nanmax(nz_cnt),1); % Create N position matrix for cell centres
    msh(ct).X=Pn(1,ct)+T(1,ct)*msh(ct).N; % Create X position matrix for cell centres
    msh(ct).Y=Pn(2,ct)+T(2,ct)*msh(ct).N; % Create Y position matrix for cell centres
    msh(ct).Z=bsxfun(@plus,bsxfun(@times,1-msh(ct).Sig,bsxfun(@minus,d_ncntr',eta_cnt)),eta_cnt); % Compute time-varying Z position for cell centers
    msh(ct).T=repmat(time_cnt,size(msh(ct).X,1),size(msh(ct).X,2)); % Create Time matrix for cells

    
    % Determine Sigma coordinates of vertices of cells
%     Maxsig=repmat(maxsig,size(msh(ct).Z,1),1);
    SLmin=cumsum([sigmin'+dsig_left'*0; repmat(dsig_left',nanmax(nz_cnt)-1,1)]); % lower left
    SLmax=cumsum([sigmin'+dsig_left'; repmat(dsig_left',nanmax(nz_cnt)-1,1)]); % upper left
    SMmin=cumsum([sigmin'+dsig_cnt'*0; repmat(dsig_cnt',nanmax(nz_cnt)-1,1)]); % lower central
    SMmax=cumsum([sigmin'+dsig_cnt'; repmat(dsig_cnt',nanmax(nz_cnt)-1,1)]); % upper central
    SRmin=cumsum([sigmin'+dsig_right'*0; repmat(dsig_right',nanmax(nz_cnt)-1,1)]); % lower right
    SRmax=cumsum([sigmin'+dsig_right'; repmat(dsig_right',nanmax(nz_cnt)-1,1)]); % upper right

    LN=repmat(nbnds(1:end-1),size(SLmin,1),1); % n of left vertices (same for upper and lower)
    MN=repmat(ncntr,size(SLmin,1),1); % n of central vertices (same for upper and lower)
    RN=repmat(nbnds(2:end),size(SLmin,1),1); % n of right vertices (same for upper and lower)
    nz_cnt(isnan(nz_cnt))=0; % use number of cells in vertical to generate fgood
    msh(ct).p.fgood=bsxfun(@le,cumsum(ones(size(msh(ct).Sig)),1),nz_cnt'); % fgood masks the data that is in the section
    msh(ct).p.N=[LN(msh(ct).p.fgood)';... % Store patch edges N coordinate
                 MN(msh(ct).p.fgood)';...
                 RN(msh(ct).p.fgood)';...
                 RN(msh(ct).p.fgood)';...
                 MN(msh(ct).p.fgood)';...
                 LN(msh(ct).p.fgood)';...
                 LN(msh(ct).p.fgood)']; % Generate matrix with n positions as is needed for patch (each column is one cell, each row a vertex)
    msh(ct).p.X=Pn(1,ct)+T(1,ct)*msh(ct).p.N; % Compute X coordinates of patch vertices
    msh(ct).p.Y=Pn(2,ct)+T(2,ct)*msh(ct).p.N; % Compute Y coordinates of patch vertices
    msh(ct).p.Sig=[SLmin(msh(ct).p.fgood)';... % Store patch edges Sigma coordinates
                 SMmin(msh(ct).p.fgood)';...
                 SRmin(msh(ct).p.fgood)';...
                 SRmax(msh(ct).p.fgood)';...
                 SMmax(msh(ct).p.fgood)';...
                 SLmax(msh(ct).p.fgood)';...
                 SLmin(msh(ct).p.fgood)']; % Generate matrix with z positions as is needed for patch (each column is one cell, each row a vertex)
     [~, nidx]=ind2sub(size(msh(ct).Sig),find(msh(ct).p.fgood)); % N index of patch edges to get bed elevation at patch corners
     ZBED= [ d_nbnds(nidx)';... % Get bed elevation for patch vertices
             d_ncntr(nidx)';...
             d_nbnds(nidx+1)';...
             d_nbnds(nidx+1)';...
             d_ncntr(nidx)';...
             d_nbnds(nidx)';...
             d_nbnds(nidx)'];                   
     msh(ct).p.Z=bsxfun(@plus,bsxfun(@times,1-msh(ct).p.Sig,bsxfun(@minus,ZBED,eta_cnt)),eta_cnt); % Compute time varying Z coordinate of patch vertices

    %% Match vertical cells for velocity data
    fnd=f_incs & fleft & inrange; % Select all data in current section, whithin the section and on the left side of cells
    dSig(fnd)=dsig_left(nIdx(fnd))+(dsig_cnt(nIdx(fnd))-dsig_left(nIdx(fnd)))/veldn*2.*(mn(fnd)-nbnds(nIdx(fnd))'); % Determine local vertical cell size (in sigma coordinates) for velocity estimates
    fnd=f_incs & fright & inrange; % Select all data in current section, whithin the section and on the right side 
    dSig(fnd)=dsig_cnt(nIdx(fnd))+(dsig_right(nIdx(fnd))-dsig_cnt(nIdx(fnd)))/veldn*2.*(mn(fnd)-ncntr(nIdx(fnd))');% Determine local vertical cell size (in sigma coordinates) for velocity estimates 

    fnde=f_incse & flefte & inrange_e;  % Select all data in current section, whithin the section and on the left side of cells (conventional processing)
    dSige(fnde)=dsig_left(nIdxe(fnde))+(dsig_cnt(nIdxe(fnde))-dsig_left(nIdxe(fnde)))/veldn*2.*(mne(fnde)-nbnds(nIdxe(fnde))');  % Determine local vertical cell size (in sigma coordinates) for velocity estimates (conventional processing)
    fnde=f_incse & frighte & inrange_e; % Select all data in current section, whithin the section and on the right side of cells (conventional processing)
    dSige(fnde)=dsig_cnt(nIdxe(fnde))+(dsig_right(nIdxe(fnde))-dsig_cnt(nIdxe(fnde)))/veldn*2.*(mne(fnde)-ncntr(nIdxe(fnde))'); % Determine local vertical cell size (in sigma coordinates) for velocity estimates (conventional processing)

    fnd=f_incs & inrange; % find velocity data belonging to current section, within the section
    sigIdx(fnd)=floor((mSig(fnd)-sigmin)./dSig(fnd))+1; % Calculate Index mapping to the right vertical cell
    
    fnde=f_incse & inrange_e; % find velocity data belonging to current section, within the section
    sigIdxe(fnde)=floor((mSige(fnde)-sigmin)./dSige(fnde))+1; % Calculate Index mapping to the right vertical cell

    %% Velocity processing
    cvel=vel(:,fcur,:);
    cTM1=TM1(:,fcur,:);
    cTM2=TM2(:,fcur,:);
    cTM3=TM3(:,fcur,:);
    % Collect all velocity data
    fnd=f_incs & sigIdx<=size(msh(ct).Z,1) & sigIdx>0; % Find all data belonging to current section and whithin a cell
    if proxim>0 % if a proximity is given
        fnd=fnd & abs(ms)<=proxim; % only include data within proxim distance from the section
    end
    siz=[size(msh(ct).Z,1) size(msh(ct).Z,2), numel(time_cnt)];
    datacol=cellfun(@(x,y,z,w) [x y z w],... % Collect velocity data and transformation matrix terms for each cell in mesh
        accumarray({sigIdx(fnd) nIdx(fnd) tIdx(fnd)},cvel(fnd),siz,@(x) {x},{}),...
        accumarray({sigIdx(fnd) nIdx(fnd) tIdx(fnd)},cTM1(fnd),siz,@(x) {x},{}),...
        accumarray({sigIdx(fnd) nIdx(fnd) tIdx(fnd)},cTM2(fnd),siz,@(x) {x},{}),...
        accumarray({sigIdx(fnd) nIdx(fnd) tIdx(fnd)},cTM3(fnd),siz,@(x) {x},{}),'UniformOutput',false); % Collect beam velocity data and corresponding transformation matrix terms
    
    % Estimate velocity (+ std, mse, Nvel)
   [msh(ct).vel, msh(ct).rsq, msh(ct).S, msh(ct).Nvel] =cellfun(@(x) estvel(x,rmresid),datacol,'UniformOutput',false); % Estimate velocity with least squares estimates (see function estvel.m for procedure)
    msh(ct).vel=squeeze(cell2mat(msh(ct).vel)); % reshape velocity data
    msh(ct).rsq=squeeze(cell2mat(msh(ct).rsq)); % reshape mean squared error data
    msh(ct).Nvel=squeeze(cell2mat(msh(ct).Nvel)); % reshape Number of valid samples data
    msh(ct).S=squeeze(cell2mat(msh(ct).S));
  
    % progressive velocity processing
    if progflag % if progressive processing is needed
       progdatacol=cell([siz nanmax(tid(ct,:))]); % Initialize variable to collect velocity and tr. matrix data
        for ccr=nanmin(tid(ct,:)):nanmax(tid(ct,:)) % loop over all crossings
           fnd=bsxfun(@and,tid(ct,fcur)<=ccr, f_incs & sigIdx<=size(msh(ct).Z,1) & sigIdx>0); % find data in current section, in all crossing up to current, and belonging to any cell
            if proxim>0 % if a proximity is given
                fnd=fnd & abs(ms)<=proxim; % only include data within proxim distance from the section
            end
            if ~any(fnd(:)), continue, end;
           progdatacol(:,:,1:siz(3),ccr)=cellfun(@(x,y,z,w) [x y z w],... %  collect velocity data and 
                accumarray({sigIdx(fnd) nIdx(fnd) tIdx(fnd)},cvel(fnd),siz,@(x) {x},{}),...
                accumarray({sigIdx(fnd) nIdx(fnd) tIdx(fnd)},cTM1(fnd),siz,@(x) {x},{}),...
                accumarray({sigIdx(fnd) nIdx(fnd) tIdx(fnd)},cTM2(fnd),siz,@(x) {x},{}),...
                accumarray({sigIdx(fnd) nIdx(fnd) tIdx(fnd)},cTM3(fnd),siz,@(x) {x},{}),'UniformOutput',false); % Collect beam velocity data and corresponding transformation matrix terms
        end
        [msh(ct).progvel, msh(ct).progrsq, msh(ct).progS, msh(ct).progNvel]=cellfun(@(x) estvel(x,rmresid), progdatacol,'UniformOutput',false); % Estimate velocity with least squares estimates (see function estvel.m for procedure)
        msh(ct).progvel=squeeze(cell2mat(msh(ct).progvel)); % reshape velocity data
        msh(ct).progrsq=squeeze(cell2mat(msh(ct).progrsq)); % reshape mean squared error data
        msh(ct).progNvel=squeeze(cell2mat(msh(ct).progNvel)); % reshape Number of valid samples data
        msh(ct).progS=squeeze(cell2mat(msh(ct).progS));
    end
%     

    % conventional velocity processing
    cvele1=vele(:,fcur,1);
    cvele2=vele(:,fcur,2);
    cvele3=vele(:,fcur,3);  
    
    fnde=f_incse & sigIdxe<=size(msh(ct).Z,1) & sigIdxe>0; % find data whithin section
 
    siz=[size(msh(ct).Z,1) size(msh(ct).Z,2), numel(time_cnt)];
    datacol=cellfun(@(x,y,z) [x y z],... % Collect velocity data and transformation matrix terms for each cell in mesh
        accumarray({sigIdxe(fnde) nIdxe(fnde) tIdx(fnde)},cvele1(fnde),siz,@(x) {x},{nan}),...
        accumarray({sigIdxe(fnde) nIdxe(fnde) tIdx(fnde)},cvele2(fnde),siz,@(x) {x},{nan}),...
        accumarray({sigIdxe(fnde) nIdxe(fnde) tIdx(fnde)},cvele3(fnde),siz,@(x) {x},{nan}),'UniformOutput',false);
    msh(ct).vele=cellfun(@(x) shiftdim(nanmean(x,1),-3),datacol,'UniformOutput',false);
    msh(ct).Se=cellfun(@(x) shiftdim(cov(x(all(isfinite(x),2),:)).*eye(3)./sum(all(isfinite(x),2)),-3),datacol,'UniformOutput',false);
    msh(ct).vele=squeeze(cell2mat(msh(ct).vele));
    msh(ct).Se=squeeze(cell2mat(msh(ct).Se));
    msh(ct).Nvele=cellfun(@(x) shiftdim(sum(all(isfinite(x),2)),-3),datacol,'UniformOutput',false);
    msh(ct).Nvele=squeeze(cell2mat(msh(ct).Nvele));
    
    if progflag % if progressive processing is needed
       progdatacol=cell([siz nanmax(tid(ct,:))]); % Initialize variable to collect velocity and tr. matrix data
        for ccr=nanmin(tid(ct,:)):nanmax(tid(ct,:)) % loop over all crossings
           fnde=bsxfun(@and,tid(ct,fcur)<=ccr, f_incse & sigIdxe<=size(msh(ct).Z,1) & sigIdxe>0); % find data in current section, in all crossing up to current, and belonging to any cell
            if proxim>0 % if a proximity is given
                fnde=fnde & abs(mse)<=proxim; % only include data within proxim distance from the section
            end
            if ~any(fnde(:)), continue, end;
           progdatacol(:,:,1:siz(3),ccr)=cellfun(@(x,y,z) [x y z],... %  collect velocity data and 
                accumarray({sigIdxe(fnde) nIdxe(fnde) tIdx(fnde)},cvele1(fnde),siz,@(x) {x},{nan}),...
                accumarray({sigIdxe(fnde) nIdxe(fnde) tIdx(fnde)},cvele2(fnde),siz,@(x) {x},{nan}),...
                accumarray({sigIdxe(fnde) nIdxe(fnde) tIdx(fnde)},cvele3(fnde),siz,@(x) {x},{nan}),'UniformOutput',false);
        end
        msh(ct).progvele=cellfun(@(x) shiftdim(nanmean(x,1),-3),progdatacol,'UniformOutput',false);
        msh(ct).progvele=squeeze(cell2mat(msh(ct).progvele));        
        msh(ct).progSe=cellfun(@(x) shiftdim(cov(x(all(isfinite(x),2),:)).*eye(3)./sum(all(isfinite(x),2)),-4),progdatacol,'UniformOutput',false);
        msh(ct).progSe=squeeze(cell2mat(msh(ct).progSe));
        msh(ct).progNvele=cellfun(@(x) shiftdim(sum(all(isfinite(x),2)),-3),progdatacol,'UniformOutput',false);
        msh(ct).progNvele=squeeze(cell2mat(msh(ct).progNvele));
    end

    
    %% Change matrix order definition
    fgood=find(msh(ct).p.fgood);
    nels=numel(msh(ct).Z);

    fgood_3=bsxfun(@plus,(0:2)*nels,fgood);
    fgood_9=bsxfun(@plus,(0:8)*nels,fgood);

    % compute new ordering
    [colm, rwm]=meshgrid(1:size(msh(ct).Z,2),1:size(msh(ct).Z,1));
    maxrow=accumarray(colm(fgood),rwm(fgood),[size(msh(ct).Z,2) 1],@max,nan)';
    rwm=bsxfun(@minus,maxrow,rwm)+1;
    fgoodr=sub2ind(size(msh(ct).Z),rwm(fgood),colm(fgood));
    fgoodr_3=bsxfun(@plus,(0:2)*nels,fgoodr);
    fgoodr_9=bsxfun(@plus,(0:8)*nels,fgoodr);
    
    % reorder data
    % init
    [tZ, tNvel, tNvele, tSig]=deal(nan(size(msh(ct).Z)));
    [tvel, tvele]=deal(nan(size(msh(ct).vel)));
    [tS, tSe]=deal(nan(size(msh(ct).S)));
    %scalar
    tZ(fgoodr)=msh(ct).Z(fgood); msh(ct).Z=tZ;
    tSig(fgoodr)=msh(ct).Sig(fgood); msh(ct).Sig=tSig;
    tNvel(fgoodr)=msh(ct).Nvel(fgood); msh(ct).Nvel=tNvel;
    tNvele(fgoodr)=msh(ct).Nvele(fgood); msh(ct).Nvele=tNvele;
    %vector
    tvel(fgoodr_3)=msh(ct).vel(fgood_3); msh(ct).vel=tvel;
    tvele(fgoodr_3)=msh(ct).vele(fgood_3); msh(ct).vele=tvele;
    %tensor
    tS(fgoodr_9)=msh(ct).S(fgood_9); msh(ct).S=tS;
    tSe(fgoodr_9)=msh(ct).Se(fgood_9); msh(ct).Se=tSe;
    
    msh(ct).p.fgood=fgoodr;
    msh(ct).p.fgood_3=fgoodr_3;
    msh(ct).p.fgood_9=fgoodr_9;
    
    if progflag
        nprog=size(msh(ct).progvel,3);
        
        mult3=reshape(0:nprog*3-1,[1 nprog 3]);
        mult9=reshape(0:nprog*9-1,[1 nprog 3 3]);      
        progfgood=bsxfun(@plus,(0:nprog-1)*nels,fgood);
        progfgood_3=bsxfun(@plus,mult3*nels,fgood);
        progfgood_9=bsxfun(@plus,mult9*nels,fgood);
        progfgoodr=bsxfun(@plus,(0:nprog-1)*nels,fgoodr);
        progfgoodr_3=bsxfun(@plus,mult3*nels,fgoodr);
        progfgoodr_9=bsxfun(@plus,mult9*nels,fgoodr);
        
        % init
        tNvel=deal(nan(size(msh(ct).progNvel)));
        tvel=deal(nan(size(msh(ct).progvel)));
        tS=deal(nan(size(msh(ct).progS)));
        tNvele=deal(nan(size(msh(ct).progNvele)));
        tvele=deal(nan(size(msh(ct).progvele)));
        tSe=deal(nan(size(msh(ct).progSe)));
        %scalar
        tNvel(progfgoodr)=msh(ct).progNvel(progfgood); msh(ct).progNvel=tNvel;
        tNvele(progfgoodr)=msh(ct).progNvele(progfgood); msh(ct).progNvele=tNvele;
        %vector
        tvel(progfgoodr_3)=msh(ct).progvel(progfgood_3); msh(ct).progvel=tvel;
        tvele(progfgoodr_3)=msh(ct).progvele(progfgood_3); msh(ct).progvele=tvele;
        %tensor
        tS(progfgoodr_9)=msh(ct).progS(progfgood_9); msh(ct).progS=tS;
        tSe(progfgoodr_9)=msh(ct).progSe(progfgood_9); msh(ct).progSe=tSe;

        msh(ct).p.progfgood=progfgoodr;
        msh(ct).p.progfgood_3=progfgoodr_3;
        msh(ct).p.progfgood_9=progfgoodr_9;
    end
    
    %% Remove data with high std
    fgood_diag=bsxfun(@plus,[0 4 8]*nels,msh(ct).p.fgood);
    if progflag
        progfgood_diag=cat(3,msh(ct).p.progfgood_9(:,:,1,1),...
                             msh(ct).p.progfgood_9(:,:,2,2),...
                             msh(ct).p.progfgood_9(:,:,3,3));
    end
    if nstd>0
		S=msh(ct).S(fgood_diag);
		Se=msh(ct).Se(fgood_diag);
        fbad=msh(ct).p.fgood(any(bsxfun(@gt,S,nstd^2*nanmedian(S,1)),2));
        fbade=msh(ct).p.fgood(any(bsxfun(@gt,Se,nstd^2*nanmedian(Se,1)),2));
        fbad_3=bsxfun(@plus,(0:2)*nels,fbad);
        fbade_3=bsxfun(@plus,(0:2)*nels,fbade);
        fbad_9=bsxfun(@plus,(0:8)*nels,fbad);
        fbade_9=bsxfun(@plus,(0:8)*nels,fbade);
        msh(ct).vel(fbad_3)=nan;
        msh(ct).S(fbad_9)=nan;
        msh(ct).vele(fbade_3)=nan;
        msh(ct).Se(fbade_9)=nan;
        if progflag
			S=msh(ct).progS(progfgood_diag);
            fbad=msh(ct).p.progfgood(any(bsxfun(@gt,S,nstd^2*S),3));
            mult3=reshape((0:2)*nprog,[1 3]);
            mult9=reshape((0:8)*nprog,[1 3 3]);      
            fbad_3=bsxfun(@plus,mult3*nels,fbad);
            fbad_9=bsxfun(@plus,mult9*nels,fbad);
            msh(ct).progvel(fbad_3)=nan;
            msh(ct).progS(fbad_9)=nan;
			Se=msh(ct).progSe(progfgood_diag);
            fbade=msh(ct).p.progfgood(any(bsxfun(@gt,S,nstd^2*Se),3));
            fbade_3=bsxfun(@plus,mult3*nels,fbade);
            fbade_9=bsxfun(@plus,mult9*nels,fbade);
            msh(ct).progvele(fbade_3)=nan;
            msh(ct).progSe(fbade_9)=nan;
			msh(ct).p.progfgood_diag=progfgood_diag;
        end
    end
    msh(ct).p.fgood_diag=fgood_diag;
    
    %% Perform rotations
    % Section averaged flow direction
    msh(ct).cs.dir=atan2(nanmean(nanmean(msh(ct).vel(:,:,2,:),1),2),nanmean(nanmean(msh(ct).vel(:,:,1,:),1),2)); % Transect averaged flow direction
    msh(ct).cs.dire=atan2(nanmean(nanmean(msh(ct).vele(:,:,2,:),1),2),nanmean(nanmean(msh(ct).vele(:,:,1,:),1),2)); % Transect averaged flow direction
    csrot=[cos(msh(ct).cs.dir) sin(msh(ct).cs.dir) 0;...
         -sin(msh(ct).cs.dir) cos(msh(ct).cs.dir) 0;...
         0 0 1];
    csrote=[cos(msh(ct).cs.dire) sin(msh(ct).cs.dire) 0;...
         -sin(msh(ct).cs.dire) cos(msh(ct).cs.dire) 0;...
         0 0 1];

    st1=matmult(msh(ct).S,shiftdim(csrot',-2),[3 4]);
    msh(ct).cs.S=matmult(shiftdim(csrot,-2),st1,[3 4]);
    msh(ct).cs.vel=matmult(shiftdim(csrot,-2),msh(ct).vel,[3 4]);
    
    st1e=matmult(msh(ct).Se,shiftdim(csrote',-2),[3 4]);
    msh(ct).cs.Se=matmult(shiftdim(csrote,-2),st1e,[3 4]);
    msh(ct).cs.vele=matmult(shiftdim(csrote,-2),msh(ct).vele,[3 4]);

    
    if progflag     
        st1=matmult(msh(ct).progS,shiftdim(csrot',-3),[4 5]);   
        msh(ct).cs.progS=matmult(shiftdim(csrot,-3),st1,[4 5]);
        msh(ct).cs.progvel=matmult(shiftdim(csrot,-3),msh(ct).progvel,[4 5]);
        st1=matmult(msh(ct).progSe,shiftdim(csrote',-3),[4 5]);   
        msh(ct).cs.progSe=matmult(shiftdim(csrote,-3),st1,[4 5]);
        msh(ct).cs.progvele=matmult(shiftdim(csrote,-3),msh(ct).progvele,[4 5]);
    end
    
    % Vertically averaged flow direction
    msh(ct).da.dir=atan2(nanmean(msh(ct).vel(:,:,2,:),1),nanmean(msh(ct).vel(:,:,1,:),1)); % depth averaged flow direction
    msh(ct).da.dire=atan2(nanmean(msh(ct).vele(:,:,2,:),1),nanmean(msh(ct).vele(:,:,1,:),1)); % depth averaged flow direction
    darot=cat(3,cat(4,cos(msh(ct).da.dir),sin(msh(ct).da.dir),zeros(size(msh(ct).da.dir))),...
                cat(4,-sin(msh(ct).da.dir),cos(msh(ct).da.dir),zeros(size(msh(ct).da.dir))),...
                cat(4,zeros(size(msh(ct).da.dir)), zeros(size(msh(ct).da.dir)),ones(size(msh(ct).da.dir))));
    darote=cat(3,cat(4,cos(msh(ct).da.dire),sin(msh(ct).da.dire),zeros(size(msh(ct).da.dir))),...
                cat(4,-sin(msh(ct).da.dire),cos(msh(ct).da.dire),zeros(size(msh(ct).da.dir))),...
                cat(4,zeros(size(msh(ct).da.dir)), zeros(size(msh(ct).da.dir)),ones(size(msh(ct).da.dir))));
    st1=matmult(msh(ct).S,permute(darot,[1 2 4 3]),[3 4]);
    msh(ct).da.S=matmult(darot,st1,[3 4]);
    msh(ct).da.vel=matmult(darot,msh(ct).vel,[3 4]);
    st1e=matmult(msh(ct).Se,permute(darote,[1 2 4 3]),[3 4]);
    msh(ct).da.Se=matmult(darote,st1e,[3 4]);
    msh(ct).da.vele=matmult(darote,msh(ct).vele,[3 4]);
    if progflag
        darot=permute(darot,[1 2 5 3 4]);
        st1=matmult(msh(ct).progS,permute(darot,[1 2 3 5 4]),[4 5]);
        msh(ct).da.progS=matmult(darot,st1,[4 5]);
        msh(ct).da.progvel=matmult(darot,msh(ct).progvel,[4 5]);
        darote=permute(darote,[1 2 5 3 4]);
        st1=matmult(msh(ct).progSe,permute(darote,[1 2 3 5 4]),[4 5]);
        msh(ct).da.progSe=matmult(darote,st1,[4 5]);
        msh(ct).da.progvele=matmult(darote,msh(ct).progvele,[4 5]);
    end
    
    % Cross-section direction
    msh(ct).sec.dir=atan2(N(2),N(1));
    secrot=[cos(msh(ct).sec.dir) sin(msh(ct).sec.dir) 0;...
         -sin(msh(ct).sec.dir) cos(msh(ct).sec.dir) 0;...
         0 0 1];
    st1=matmult(msh(ct).S,shiftdim(secrot',-2),[3 4]);
    msh(ct).sec.S=matmult(shiftdim(secrot,-2),st1,[3 4]);
    msh(ct).sec.vel=matmult(shiftdim(secrot,-2),msh(ct).vel,[3 4]);
    st1e=matmult(msh(ct).Se,shiftdim(secrot',-2),[3 4]);
    msh(ct).sec.Se=matmult(shiftdim(secrot,-2),st1e,[3 4]);
    msh(ct).sec.vele=matmult(shiftdim(secrot,-2),msh(ct).vele,[3 4]);   
    if progflag     
        st1=matmult(msh(ct).progS,shiftdim(secrot',-3),[4 5]);   
        msh(ct).sec.progS=matmult(shiftdim(secrot,-3),st1,[4 5]);
        msh(ct).sec.progvel=matmult(shiftdim(secrot,-3),msh(ct).progvel,[4 5]);
        st1=matmult(msh(ct).progSe,shiftdim(secrot',-3),[4 5]);   
        msh(ct).sec.progSe=matmult(shiftdim(secrot,-3),st1,[4 5]);
        msh(ct).sec.progvele=matmult(shiftdim(secrot,-3),msh(ct).progvele,[4 5]);
    end

    
end



function [vel, rsq, S, nin]=estvel(in,rmresid)
% estimates cartesian velocities given beam velocities and correspoding
% transformation matrix terms. IN is Nx4, N being the amount of beam
% velocity samples. Each row is composed of a beam velocity and the three
% terms from the transformation matrix

in(any(isnan(in),2),:)=[]; % Throw out bad data
rmvel=true; % initialize variable holding indices to outliers

while any(rmvel) % iterate until no outliers are found
    if isempty(in) || rank(in(:,2:4))<3 % if no velocity is available or matrix is rank deficient
        vel=nan(3,1); rsq=nan; S=nan(3,3); % return nans
        break; % Exit loop
    else
        [vel, ~, ~, S]=lscov(in(:,2:4),in(:,1)); % Perform least squares estimate of velocities
        res=abs(in(:,1)-in(:,2:4)*vel); % Determine residuals in beam velocity
        rsq=1-sum(res.^2)./sum((in(:,1)-mean(in(:,1))).^2);
    end
    if rmresid==0, break, end % if rmvel == 0, no outlier removal needed, exit loop
    rmvel=res>rmresid*median(res); % Detect outliers from residuals
    in(rmvel,:)=[]; % Remove outliers
end
nin=shiftdim(size(in,1),-3); % compute amount of points used for final velocity computation
vel=shiftdim(vel,-4); % shift dimension to output 1x1x1x1x3 vector
rsq=shiftdim(rsq,-4); % shift dimension to output 1x1x1x1x3 vector
S=shiftdim(S,-4); % shift dimension to output 1x1x1x1x3x3 vector



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
