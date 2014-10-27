function msh=procTrans2(adcp,tid,varargin)
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


%    Copyright 2013,2014 Bart Vermeulen
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

% TODO:
%   Add possibility to define maximum_z or maximum_sigma for mesh
%   Smooth bathymetry interpolator


%% Handle input
assert(isADCPstruct(adcp)); % Check first input is an ADCP structure
assert(size(tid,2)==size(adcp.VEL,2)); % Ensure size of second dimension for the transect id's matches the number of ensembles in the data

nens=size(adcp.VEL,2);

P=inputParser;
P.FunctionName='procTrans';
P.addParamValue('DepthTransducer',0.3,@(x) isscalar(x) && isnumeric(x) && x>0) % depth of transducer
P.addParamValue('DeltaN',5,@(x) isscalar(x) && isnumeric(x) && x>0) % N resolution of mesh
P.addParamValue('DeltaZ',1,@(x) isscalar(x) && isnumeric(x) && x>0) % target Z resolution of mesh
P.addParamValue('Eta',0,@(x) isvector(x) && isnumeric(x) && (numel(x)==1 || numel(x)==nens)) % eta (water level)
P.addParamValue('MinimumSigma',0.06,@(x) isscalar(x) && isnumeric(x) && x>=0 && x <=1) % Minimum sigma for meshing (to account for side lobes)
P.addParamValue('CumulateCrossings',false,@(x) islogical(x) && isscalar(x)) % Set whether to lump consecutive crossings (only works if crossings are indicated in tid)
P.addParamValue('ConventionalProcessing',false,@(x) islogical(x) && isscalar(x)) % Set whether to use conventional processing
P.addParamValue('TopMeshLowestEta',false,@(x) islogical(x) && isscalar(x)) % Set whether top of mesh should be set at lowest water level (or top of data if this is lower than lowest water level)
P.addParamValue('ConstantZetaMesh',false,@(x) islogical(x) && isscalar(x)) % Set whether to keep sigma constant. If false z is kept constant
P.addParamValue('RemoveOutliers',0,@(x) isscalar(x) && isnumeric(x) && x>=0); % Remove outliers when their residuals exceed this amount of times the median residual in velocity inversion
P.addParamValue('Proximity',0,@(x) isscalar(x) && isnumeric(x) && x>=0); % Set the maximum distance from the cross-section for data to be included in the calculation
P.addParamValue('StdFiltering',6,@(x) isscalar(x) && isnumeric(x) && x>=0); % Remove velocity with a standard deviation of the estimate exceeding this amount of times the median standard deviation
P.addParamValue('ShipReference','bt',@(x) ischar(x) && any(strcmpi(x,{'bt','gps','gpsbt','btgps'}))); % Select which boat velocity calculation to use
P.addParamValue('Pusr',[], @(x) isnumeric(x)); % User set the central position of a transect (can be usefull when data density is not uniform on transect
P.addParamValue('ModelU_t',1,@(x) isnumeric(x) && ~isempty(x) && (isscalar(x) || (ismatrix(x) && size(x,2)==nens))); % Velocity model for u
P.addParamValue('ModelV_t',1,@(x) isnumeric(x) && ~isempty(x) && (isscalar(x) || (ismatrix(x) && size(x,2)==nens))); % Velocity model for v
P.addParamValue('ModelW_t',1,@(x) isnumeric(x) && ~isempty(x) && (isscalar(x) || (ismatrix(x) && size(x,2)==nens))); % Velocity model for w
P.addParamValue('RotatePars',true,@(x) isscalar(x) && islogical(x));
P.addParamValue('EnableDebugging',false,@(x) isscalar(x) && islogical(x));

P.parse(varargin{:}); % parse the input


% Retrieve all input
depthtransd=P.Results.DepthTransducer;
veldn=P.Results.DeltaN;
veldz=P.Results.DeltaZ;
eta=P.Results.Eta;
sigmin=P.Results.MinimumSigma;
rmresid=P.Results.RemoveOutliers;
proxim=P.Results.Proximity;
nstd=P.Results.StdFiltering;
maxz_mineta=P.Results.TopMeshLowestEta;
constant_sigma_mesh=~P.Results.ConstantZetaMesh;
time=datenum(adcp.timeV)'; % Get time
eta=eta(:)';    % Vectorize eta
eta=eta.*ones(1,size(numel(time),1)); % Ensure eta has same number of elements as the number of ensembles
shref=P.Results.ShipReference;
Pusr=P.Results.Pusr;
mod_ut=P.Results.ModelU_t;
mod_ut=permute(shiftdim(bsxfun(@times,mod_ut,ones(1,nens)),-2),[1 4 2 3]); % ensure size(mod_ut) is [1 nens 1 npars]
mod_vt=P.Results.ModelV_t;
mod_vt=permute(shiftdim(bsxfun(@times,mod_vt,ones(1,nens)),-2),[1 4 2 3]); % ensure size(mod_ut) is [1 nens 1 npars]
mod_wt=P.Results.ModelW_t;
mod_wt=permute(shiftdim(bsxfun(@times,mod_wt,ones(1,nens)),-2),[1 4 2 3]); % ensure size(mod_ut) is [1 nens 1 npars]
npars_u=size(mod_ut,4);
npars_v=size(mod_vt,4);
npars_w=size(mod_wt,4);
npars=npars_u+npars_v+npars_w;
IsCumulative=P.Results.CumulateCrossings;
IsConventional=P.Results.ConventionalProcessing;
rotate=P.Results.RotatePars;
fdebug=P.Results.EnableDebugging;

%% Get ADCP positioning data
[xadcp,yadcp]=utmADCP(adcp); % Get adcp position
misal=getExtMisalign(adcp); % Get beam 3 misalignment

% depth location
[xb,yb,zb]=depthADCP(adcp,'Beam3Misalign',misal); % Get offsets from ADCP to the bed
zb=bsxfun(@plus,zb-depthtransd,eta); % transform vertical offset to bed level with respect to eta=0
xb=bsxfun(@plus,xadcp',xb); % transform offsets to UTM coordinate system for horizontal
yb=bsxfun(@plus,yadcp',yb); % transform offsets to UTM coordinate system for horizontal
if IsConventional
    xb=nanmean(xb,3); % determine beam average for conventional velocity processing
    yb=nanmean(yb,3); % determine beam average for conventional velocity processing
    zb=nanmean(zb,3); % determine beam average for conventional velocity processing
end


% velocity location
[xvel,yvel,zvel]=mapADCP(adcp,'IsUpward', false, 'Beam3Misalign',misal); % get offsets from adcp to velocity locations
zvel=bsxfun(@plus,zvel-depthtransd,eta); % transform vertical offset to z position of velocity with respect to eta=0
xvel=bsxfun(@plus,xadcp',xvel); % transform offsets to vel data to UTM coordinate system for horizontal
yvel=bsxfun(@plus,yadcp',yvel); % transform offsets to vel data to UTM coordinate system for horizontal
if IsConventional
    xvel=nanmean(xvel,3); % determine beam average for conventional velocity processing
    yvel=nanmean(yvel,3); % determine beam average for conventional velocity processing
    zvel=nanmean(zvel,3); % determine beam average for conventional velocity processing
end


% velocity data
[adcp.VEL,adcp.btvel] = filterADCP(adcp,'','filterBT',true); % Filter velocity
if IsConventional
    gpsvel=getGPSvel(adcp,'e');
    [vel,btvel]=corADCP(adcp,'e','UseExtHeading',true,'Beam3Misalign',misal); % transform to earth velocity
else
    gpsvel=getGPSvel(adcp,'b'); 
    [vel, btvel]=corADCP(adcp,'b','UseExtHeading',true,'Beam3Misalign',misal); % transform to beam velocity
    if fdebug
        adcp2=adcp;
        [adcp2.VEL,adcp2.btvel]=corADCP(adcp2,'e','UseExtHeading',true,'Beam3Misalign',misal); % transform to earth velocity;
        [bvel,bbtvel]=corADCP(adcp2,'b','UseExtHeading',true,'Beam3Misalign',misal,'ForceOrigin','e'); % transform to earth velocity;
        imagesc(adcp.VEL(:,:,1)-bvel(:,:,1))
        colorbar
    end
end

fbad_bt=any(isnan(btvel),2);
zb(:,fbad_bt,:)=nan;

if strcmpi(shref,'gps')
    btvel=gpsvel;
elseif strcmpi(shref,'btgps')
    btvel(fbad_bt,:)=gpsvel(fbad_bt,:);
elseif strcmpi(shref,'gpsbt')
    fbad_gps=any(isnan(gpsvel),2);
    gpsvel(fbad_gps,:)=btvel(fbad_gps,:);
    btvel=gpsvel;
end

vel=vel-repmat(shiftdim(btvel,-1),[size(vel,1),1,1]); % remove bt from normal vel

if ~IsConventional 
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
else
    [TM1, TM2, TM3]=deal(zeros(size(vel)));
    TM1(:,:,1)=1;
    TM2(:,:,2)=1;
    TM3(:,:,3)=1;
end


%% Project positions to section
% Initialize
nsec=size(tid,1);
nens=size(vel,2);
msh=repmat(struct(),[nsec,1]); % Initialize final structure for output
T=nan(2,nsec); % Unit vector tangential to section
N=nan(2,nsec); % Unit vector orthogonal to section 
Pm=nan(2,nsec); % Average section location vector
Pn=nan(2,nsec); % Average section location vector


%% create time index (mapping to right time cell)
eta = eta .* ones(1,nens); % make sure eta is row vector with same number of elements as number of ensembles

% Start processing
for ct=1:size(tid,1) % For all sections  
    fcur=tid(ct,:)>0; % mask, indicating data belonging to current section
    nrepeat=nanmax(tid(ct,:));
    % compute eta and time for each crossing (or cumulated crossings)
    msh(ct).eta=accumarray(tid(ct,fcur)',eta(fcur)',[],@nanmean,nan); % Compute average eta for each crossing
    msh(ct).maxeta=accumarray(tid(ct,fcur)',eta(fcur)',[],@nanmax,nan); % Compute maximum eta for each crossing
    msh(ct).mineta=accumarray(tid(ct,fcur)',eta(fcur)',[],@nanmin,nan); % Compute minimum eta for each crossing
    msh(ct).time=accumarray(tid(ct,fcur)',time(fcur)',[],@nanmean,nan); % Compute average eta for each crossing
    if IsCumulative
        maxprev=-inf;
        minprev=inf;
        for cs=1:nsec
            msh(ct).maxeta(cs)=max(msh(ct).maxeta(cs),maxprev);
            maxprev=msh(ct).maxeta(cs);
            msh(ct).mineta(cs)=min(msh(ct).mineta(cs),minprev);
            minprev=msh(ct).mineta(cs);
        end
        msh(ct).eta=cumsum(msh(ct).eta)./(1:numel(msh(ct).eta))';% Compute cumulative average if processing is cumulative
    end
    
    % Compute relevant vectors
    P=[xadcp(fcur)';yadcp(fcur)']; % ADCP position vectors
    if (~isempty(Pusr))
        Pm = mean(Pusr(:,:,ct),2);
        T(:,ct) = diff(Pusr(:,:,ct),[],2);
        T(:,ct) = T(:,ct)/norm(T(:,ct));
    else
        Pm(:,ct)=nanmean(P,2); % Average ADCP position vector
        T(:,ct)=princdir(P); % Section tangential vector, determined as largest eigenvector of P
    end
    N(:,ct)=[T(2,ct); -T(1,ct)]; % Section orthogonal vector (orthogonal to T)
    
    % Project positions (inner product with unit vectors)
    nb=T(1,ct)*(xb(:,fcur,:)-Pm(1,ct))+T(2,ct)*(yb(:,fcur,:)-Pm(2,ct)); % depth position, n component
    sb=N(1,ct)*(xb(:,fcur,:)-Pm(1,ct))+N(2,ct)*(yb(:,fcur,:)-Pm(2,ct)); % depth position, s component
    nvel=T(1,ct)*(xvel(:,fcur,:)-Pm(1,ct))+T(2,ct)*(yvel(:,fcur,:)-Pm(2,ct)); % velocity position, n component
    svel=N(1,ct)*(xvel(:,fcur,:)-Pm(1,ct))+N(2,ct)*(yvel(:,fcur,:)-Pm(2,ct)); % velocity position, s component
    
    % Store important vectors
    msh(ct).Tvec=T(:,ct); % tangential vector
    msh(ct).Nvec=N(:,ct); % orthogonal vector
    msh(ct).Pm=Pm(:,ct); % mean position vector
    
    %% Calculate sigma for velocity locations
    czb=zb(:,fcur,:);
    
    % Finite depth and velocity masks
    dfg=isfinite(nb) & isfinite(sb) & isfinite(czb); % Mask for finite depth locations
    velfg=isfinite(nvel) & isfinite(svel); % Mask for finite velocity locations

    % Initialize variables
    zbvel=nan(size(nvel)); % Bed elevation at velocity locations

    % Compute Bed elevations at velocity sample locations (use accumarray instead?)
    if verLessThan('matlab','7.10')
        zbvel(velfg) = griddata(reshape(sb(dfg),[],1),reshape(nb(dfg),[],1),reshape(czb(dfg),[],1),svel(velfg),nvel(velfg),'natural'); %#ok<GRIDD>
    elseif verLessThan('matlab','8.1')
        Int=TriScatteredInterp(reshape(sb(dfg),[],1),reshape(nb(dfg),[],1),reshape(czb(dfg),[],1),'natural'); %#ok<REMFF1> % Create interpolant (natural interpolation with linear extrapolation)
        zbvel(velfg)=Int(svel(velfg),nvel(velfg)); % interpolate bed elevation at velocity locations
    else
        Int=scatteredInterpolant(reshape(sb(dfg),[],1),reshape(nb(dfg),[],1),reshape(czb(dfg),[],1),'natural','none'); % Create interpolant (natural interpolation with linear extrapolation)
        zbvel(velfg)=Int(svel(velfg),nvel(velfg)); % interpolate bed elevation at velocity locations
    end  

    % Compute sigma
    czvel=zvel(:,fcur,:);
    if numel(eta)>1, ceta=eta(fcur); else ceta=eta; end
    sigVel=nan(size(czvel));
    f_incs=bsxfun(@gt,czvel,zbvel);
    dVel=bsxfun(@minus,ceta,zbvel);
    sigVel(f_incs)=(czvel(f_incs)-zbvel(f_incs))./dVel(f_incs); % Compute sigma for velocity locations
    f_incs=sigVel<=1 & sigVel>=sigmin; % Create mask for velocity locations in the cross-section
    
    
    %% Meshing
    % Compute minimum and maximum n
    if (~isempty(Pusr))
        nmin = min(T(:,ct)'*(Pusr(:,:,ct) - Pm(:,ct)*[1 1]));
    else
    nmin=nanmin(nvel(f_incs)); % determine minimum n
    end
    Pn(:,ct)=Pm(:,ct)+T(:,ct)*nmin; % compute vector pointing to minimum n location (projected on section)
    nvel=nvel-nmin; % remove minimum n from velocity n coordinates (now they range between 0 and maximum n)
    nb=nb-nmin; % remove minimum n from depth n coordinates (now they range between 0 and maximum n)
    if (~isempty(Pusr))
        nmax = max(T(:,ct)'*(Pusr(:,:,ct) - Pn(:,ct)*[1 1]));
    else
    nmax=nanmax(nvel(f_incs)); % compute maximum n
    end
    ncells=ceil(nmax./veldn); % determine number of vertical in section

    % Preliminary computations
    nIdxd=floor((nb+veldn/4)/(veldn/2))+1; % Index to determine depth at cell boundaries and centres
    fgood=isfinite(nIdxd) & nIdxd>0 & nIdxd<=ncells*2+1; % Select data to include in depth calculation at cell edges and center
    fleft=mod(floor(nvel/veldn*2),2)==0; % Selects data on the left half of a cell
    fright=~fleft; % Selects data on the left half of a cell

    % determine n coordinate and bed_elevation of cell edges and centres (for crossing with maximum eta)
    zb_tmp=accumarray(reshape(nIdxd(fgood),[],1),reshape(czb(fgood),[],1),[ncells*2+1,1],@nanmean,nan); % Determine elevation (z) at cell centres and boundaries (where available)
    fgood_zb=isfinite(zb_tmp); % find all available elevations at cell centres and boundaries
    Int=griddedInterpolant(find(fgood_zb),zb_tmp(fgood_zb)); % create interpolant to calculate elevation where it's missing (probably at the edges)
    zb_tmp(~fgood_zb)=Int(find(~fgood_zb)); %#ok<FNDSB> % interpolate at cell centres or boundaries with missing depth

    % Compute maximum z verticals
    maxEta=max(msh(ct).maxeta); % Get maximum crossing avareged eta   
    if constant_sigma_mesh
        maxsig=nanmax(sigVel(f_incs)); % maximum Sigma measured in all data belonging to current section
        max_zb=nanmax(zb_tmp);
        maxz=max_zb+maxsig*(maxEta-max_zb); % Transform maximum sigma to maximum mesh elevation at section location
    else
        maxz=nanmax(czvel(:));
        if ~constant_sigma_mesh && maxz_mineta
            maxz=min(maxz,nanmin(msh(ct).mineta));
        end
    end
    
    % check for degenerate cells (left or right boundary exceeds maximum z)
    % and fix them (This should only be due to extrapolation and should not
    % occur anywhere except at the edges)
    ntmp=(0:0.5:ncells)*veldn; 
    if zb_tmp(1)>maxz
        ntmp(1)=ntmp(2)-abs((zb_tmp(2)-maxz)./(zb_tmp(2)-zb_tmp(1))*(veldn/2));
        zb_tmp(1)=maxz;
    end
    if zb_tmp(end)>maxz
        ntmp(end)=ntmp(end-1)+abs((zb_tmp(end-1)-maxz)./(zb_tmp(end)-zb_tmp(end-1))*(veldn/2));
        zb_tmp(end)=maxz;
    end
    if(any(zb_tmp(2:end-1)>maxz))
        warning('procTrans2:DegenerateCells','Degenerate cells found (round 1)!')
    end
    
    % compute minimum z
    d_tmp=maxEta-zb_tmp; % depth at cell boundaries and centres
    minz_tmp=zb_tmp+sigmin*d_tmp; % minimum bed level for mesh (centres and boundaries)
    if minz_tmp(1)>maxz
        ntmp(1)=ntmp(2)-abs((minz_tmp(2)-maxz)./(minz_tmp(2)-minz_tmp(1))*(veldn/2));
        minz_tmp(1)=maxz;
        zb_tmp(1)=(minz_tmp(1)-sigmin*maxEta)/(1-sigmin);
        d_tmp(1)=maxEta-zb_tmp(1);
    end
    if minz_tmp(end)>maxz
        ntmp(end)=ntmp(end-1)+abs((minz_tmp(end-1)-maxz)./(minz_tmp(end)-minz_tmp(end-1))*(veldn/2));
        minz_tmp(end)=maxz;
        zb_tmp(end)=(minz_tmp(end)-sigmin*maxEta)/(1-sigmin);
        d_tmp(end)=maxEta-zb_tmp(end);
    end
    if(any(minz_tmp(2:end-1)>maxz))
        warning('procTrans2:DegenerateCells','Degenerate cells found (round 2)!')
    end

    
    % store some stuff usefull for plotting
    msh(ct).p.nbed=ntmp; % n coordinate of both cell centres and boundaries
    msh(ct).p.xbed=Pn(1,ct)+T(1,ct)*msh(ct).p.nbed; % x coordinate of both cell centres and boundaries
    msh(ct).p.ybed=Pn(2,ct)+T(2,ct)*msh(ct).p.nbed; % y coordinate of both cell centres and boundaries
    msh(ct).p.zbed=zb_tmp'; % store the z coordinate of the bed for cell centres and boundaries

    
    d_ncntr=d_tmp(2:2:end); % get elevation at cell centres
    d_nbnds=d_tmp(1:2:end); % get average bed elevation at cell boundaries
    d_nleft=d_nbnds(1:end-1); % get average bed elevation at cell left boundary
    d_nright=d_nbnds(2:end); % get average bed elevation at cell right boundary
    
    zb_ncntr=zb_tmp(2:2:end); % get elevation at cell centres
    zb_nbnds=zb_tmp(1:2:end); % get average bed elevation at cell boundaries
    zb_nleft=zb_nbnds(1:end-1); % get average bed elevation at cell left boundary
    zb_nright=zb_nbnds(2:end); % get average bed elevation at cell right boundary
    
    minz_cnt=minz_tmp(2:2:end); % minimum z for each vertical at cell centres
    minz_bnd=minz_tmp(1:2:end); % minimum z for each vertical at cell boundaries   
    minz_left=minz_bnd(1:end-1); % minimum z for each vertical at cell left boundaries
    minz_right=minz_bnd(2:end); % minimum z for each vertical at cell right boundaries

    n_cntr=ntmp(2:2:end);
    n_bnds=ntmp(1:2:end);
    n_left=n_bnds(1:end-1);
    n_right=n_bnds(2:end);
    
    % Compute n indices
    colIdx=floor((nvel-n_right(1))./veldn)+2; % Compute N index, mapping velocity data to correct vertical

    % select points within transect range
    inrange = colIdx >0 & colIdx <= ncells;
    colIdx(colIdx < 1 & colIdx > ncells) = NaN;    
    
    
    % Compute number of cells in each vertical and cell-size (in sigma and zed)
    nz_cnt=round((maxz-minz_cnt)/veldz); % best guess number of z values in vertical giving a dz as close as possible to the given one
    dzed_cnt=(maxz-minz_cnt)./nz_cnt; % compute dz at cell centers
    dzed_left=(maxz-minz_left)./nz_cnt; % compute dz at left cell boundary
    dzed_right=(maxz-minz_right)./nz_cnt; % compute dz at right cell boundary
    
    % generate N coordinates of the mesh
    NL=repmat(n_left,nanmax(nz_cnt),1); % n of left vertices (same for upper and lower)
    NM=repmat(n_cntr,nanmax(nz_cnt),1); % n of central vertices (same for upper and lower)
    NR=repmat(n_right,nanmax(nz_cnt),1); % n of right vertices (same for upper and lower)

    % make fgood
    nz_cnt(isnan(nz_cnt))=0; % use number of cells in vertical to generate fgood
    msh(ct).p.fgood=bsxfun(@le,cumsum(ones(size(NL)),1),nz_cnt'); % fgood masks the data that is in the section  
    n_rep_trans=nanmax(tid(ct,:));
    mult_vec=reshape(0:n_rep_trans*npars-1,[1 n_rep_trans npars]);
    mult_tens=reshape(0:n_rep_trans*npars^2-1,[1 n_rep_trans npars npars]);  
    nels=numel(NL);
    fgood=find(msh(ct).p.fgood);
    msh(ct).p.progfgood=bsxfun(@plus,(0:n_rep_trans-1)*nels,fgood);
    msh(ct).p.progfgood_vec=bsxfun(@plus,mult_vec*nels,fgood);
    msh(ct).p.progfgood_tens=bsxfun(@plus,mult_tens*nels,fgood);
    col_fgood=cell(npars,1);
    for cpars=1:npars
        col_fgood{cpars}=msh(ct).p.progfgood_tens(:,:,cpars,cpars);
    end
    msh(ct).p.progfgood_tensdiag=cat(3,col_fgood{:});


    

    % generate N, X and Y coordinates for centres and patch vertices
    msh(ct).N=repmat(n_cntr,nanmax(nz_cnt),1); % Create N position matrix for cell centres
    msh(ct).X=Pn(1,ct)+T(1,ct)*msh(ct).N; % Create X position matrix for cell centres
    msh(ct).Y=Pn(2,ct)+T(2,ct)*msh(ct).N; % Create Y position matrix for cell centres
    msh(ct).p.N=[NL(msh(ct).p.fgood)';... % Store patch edges N coordinate
                 NM(msh(ct).p.fgood)';...
                 NR(msh(ct).p.fgood)';...
                 NR(msh(ct).p.fgood)';...
                 NM(msh(ct).p.fgood)';...
                 NL(msh(ct).p.fgood)';...
                 NL(msh(ct).p.fgood)']; % Generate matrix with n positions as is needed for patch (each column is one cell, each row a vertex)
    msh(ct).p.X=Pn(1,ct)+T(1,ct)*msh(ct).p.N; % Compute X coordinates of patch vertices
    msh(ct).p.Y=Pn(2,ct)+T(2,ct)*msh(ct).p.N; % Compute Y coordinates of patch vertices

    % Bed elevations at patch vertices
    [~, idx_pedges]=ind2sub(size(msh(ct).N),find(msh(ct).p.fgood)); % N index of patch edges to get bed elevation at patch corners
    ZBED= [ zb_nleft(idx_pedges)';... % Get bed elevation for patch vertices
            zb_ncntr(idx_pedges)';...
            zb_nright(idx_pedges)';...
            zb_nright(idx_pedges)';...
            zb_ncntr(idx_pedges)';...
            zb_nleft(idx_pedges)';...
            zb_nleft(idx_pedges)'];                   

    % Initialize row idx
    rowIdx=nan(size(nvel)); % Sigma index (mapping to right depth cell)

    % generate Sigma and Z coordinates of mesh   
    if constant_sigma_mesh      
        dsig_cnt=dzed_cnt./d_ncntr; % compute dsigma at cell centers
        dsig_left=dzed_left./d_nleft; % compute dsigma at left cell boundary
        dsig_right=dzed_right./d_nright; % comput dsigma at right cell boundary

         % Determine Sigma coordinates of mesh
        Scnt=cumsum([sigmin+dsig_cnt'/2; repmat(dsig_cnt',nanmax(nz_cnt)-1,1)]); % generate sigma positions of cell centres
        SLmin=cumsum([sigmin'+dsig_left'*0; repmat(dsig_left',nanmax(nz_cnt)-1,1)]); % lower left
        SLmax=cumsum([sigmin'+dsig_left'; repmat(dsig_left',nanmax(nz_cnt)-1,1)]); % upper left
        SMmin=cumsum([sigmin'+dsig_cnt'*0; repmat(dsig_cnt',nanmax(nz_cnt)-1,1)]); % lower central
        SMmax=cumsum([sigmin'+dsig_cnt'; repmat(dsig_cnt',nanmax(nz_cnt)-1,1)]); % upper central
        SRmin=cumsum([sigmin'+dsig_right'*0; repmat(dsig_right',nanmax(nz_cnt)-1,1)]); % lower right
        SRmax=cumsum([sigmin'+dsig_right'; repmat(dsig_right',nanmax(nz_cnt)-1,1)]); % upper right
    
        % reorder matrix (small sigma at the bottom)
        fgood=find(msh(ct).p.fgood);
        [colm, rwm]=meshgrid(1:size(Scnt,2),1:size(Scnt,1));
        maxrow=accumarray(colm(fgood),rwm(fgood),[size(Scnt,2) 1],@max,nan)';
        rwm=bsxfun(@minus,maxrow,rwm)+1;
        fgoodr=sub2ind(size(Scnt),rwm(fgood),colm(fgood));
        Scnt(fgoodr)=Scnt(fgood);
        SLmin(fgoodr)=SLmin(fgood);
        SLmax(fgoodr)=SLmax(fgood);
        SMmin(fgoodr)=SMmin(fgood);
        SMmax(fgoodr)=SMmax(fgood);
        SRmin(fgoodr)=SRmin(fgood);
        SRmax(fgoodr)=SRmax(fgood);
        
        % make mesh in sigma
        msh(ct).p.Sig=repmat([SLmin(msh(ct).p.fgood)';... % Store patch edges Sigma coordinates
                              SMmin(msh(ct).p.fgood)';...
                              SRmin(msh(ct).p.fgood)';...
                              SRmax(msh(ct).p.fgood)';...
                              SMmax(msh(ct).p.fgood)';...
                              SLmax(msh(ct).p.fgood)';...
                              SLmin(msh(ct).p.fgood)'],1,1,nanmax(tid(ct,:))); % Generate matrix with z positions as is needed for patch (each column is one cell, each row a vertex)
        msh(ct).Sig=Scnt;
        msh(ct).Sig(~msh(ct).p.fgood)=nan;
        msh(ct).Sig=repmat(msh(ct).Sig,1,1,nanmax(tid(ct,:)));
        
        % compute mesh in Z
        msh(ct).Z=bsxfun(@plus,bsxfun(@times,msh(ct).Sig,bsxfun(@minus,reshape(msh(ct).eta,1,1,[]),zb_ncntr')),zb_ncntr'); % z= zb+sig(eta-zb)
        msh(ct).p.Z=bsxfun(@plus,bsxfun(@times,msh(ct).p.Sig,bsxfun(@minus,reshape(msh(ct).eta,1,1,[]),ZBED)),ZBED); % z= zb+sig(eta-zb)
        
        % Make row indices for velocity
        vel_maxsig=nan(size(sigVel));
        fnd=fleft & inrange;
        vel_maxsig(fnd)=SLmax(1,colIdx(fnd))+(SMmax(1,colIdx(fnd))-SLmax(1,colIdx(fnd)))./(n_cntr(colIdx(fnd))-n_left(colIdx(fnd))).*(nvel(fnd)-n_left(colIdx(fnd))')';
        fnd=fright & inrange;
        vel_maxsig(fnd)=SMmax(1,colIdx(fnd))+(SRmax(1,colIdx(fnd))-SMmax(1,colIdx(fnd)))./(n_right(colIdx(fnd))-n_cntr(colIdx(fnd))).*(nvel(fnd)-n_cntr(colIdx(fnd))')';
        f_incs=sigVel<=vel_maxsig & sigVel >= sigmin;
        dSig=nan(size(nvel)); % Step in sigma
        fnd=f_incs & fleft & inrange; % Select all data in current section, whithin the section and on the left side of cells
        dSig(fnd)=dsig_left(colIdx(fnd))+(dsig_cnt(colIdx(fnd))-dsig_left(colIdx(fnd)))./(n_cntr(colIdx(fnd))-n_left(colIdx(fnd)))'.*(nvel(fnd)-n_left(colIdx(fnd))'); % Determine local vertical cell size (in sigma coordinates) for velocity estimates
        fnd=f_incs & fright & inrange; % Select all data in current section, whithin the section and on the right side 
        dSig(fnd)=dsig_cnt(colIdx(fnd))+(dsig_right(colIdx(fnd))-dsig_cnt(colIdx(fnd)))./(n_right(colIdx(fnd))-n_cntr(colIdx(fnd)))'.*(nvel(fnd)-n_cntr(colIdx(fnd))');% Determine local vertical cell size (in sigma coordinates) for velocity estimates 
        fnd=f_incs & inrange; % find velocity data belonging to current section, within the section
        rowIdx(fnd)=floor((sigVel(fnd)-sigmin)./dSig(fnd))+1; % Calculate Index mapping to the right vertical cell
        rowIdx(fnd)=maxrow(colIdx(fnd))-rowIdx(fnd)'+1; % Reverse indices to have data measured on top in the top of the matrix
    else
        % Determine z-coordinates of mesh
        Zcnt=cumsum([maxz-dzed_cnt'/2; repmat(-dzed_cnt',nanmax(nz_cnt)-1,1)],1);
        ZLmax=cumsum([maxz-dzed_left'*0; repmat(-dzed_left',nanmax(nz_cnt)-1,1)],1);
        ZLmin=cumsum([maxz-dzed_left'; repmat(-dzed_left',nanmax(nz_cnt)-1,1)],1);
        ZMmax=cumsum([maxz-dzed_cnt'*0; repmat(-dzed_cnt',nanmax(nz_cnt)-1,1)],1);
        ZMmin=cumsum([maxz-dzed_cnt'; repmat(-dzed_cnt',nanmax(nz_cnt)-1,1)],1);
        ZRmax=cumsum([maxz-dzed_right'*0; repmat(-dzed_right',nanmax(nz_cnt)-1,1)],1);
        ZRmin=cumsum([maxz-dzed_right'; repmat(-dzed_right',nanmax(nz_cnt)-1,1)],1);
        
        % make mesh in Z
        msh(ct).p.Z=repmat([ZLmin(msh(ct).p.fgood)';... % Store patch edges Z coordinates
                            ZMmin(msh(ct).p.fgood)';...
                            ZRmin(msh(ct).p.fgood)';...
                            ZRmax(msh(ct).p.fgood)';...
                            ZMmax(msh(ct).p.fgood)';...
                            ZLmax(msh(ct).p.fgood)';...
                            ZLmin(msh(ct).p.fgood)'],1,1,nanmax(tid(ct,:))); % Generate matrix with z positions as is needed for patch (each column is one cell, each row a vertex)
        msh(ct).Z=Zcnt;
        msh(ct).Z(~msh(ct).p.fgood)=nan;
        msh(ct).Z=repmat(msh(ct).Z,1,1,nanmax(tid(ct,:)));
        
        % Compute mesh in Sigma
        msh(ct).Sig=bsxfun(@rdivide,bsxfun(@minus,msh(ct).Z,zb_ncntr'),bsxfun(@minus,reshape(msh(ct).eta,1,1,[]),zb_ncntr')); % sig=(z-zb)/(eta-zb)
        msh(ct).p.Sig=bsxfun(@rdivide,bsxfun(@minus,msh(ct).p.Z,ZBED),bsxfun(@minus,reshape(msh(ct).eta,1,1,[]),ZBED)); % sig=(z-zb)/(eta-zb)
        
        % Make row indices for velocity
        vel_minz=nan(size(czvel));
        fnd=fleft & inrange;
        frw=sum(msh(ct).p.fgood,1);
        idx_cell=sub2ind(size(msh(ct).Z),max(1,frw(colIdx(fnd))),colIdx(fnd)'); % max is quick and dirty solution?
        vel_minz(fnd)=ZLmin(idx_cell)+(ZMmin(idx_cell)-ZLmin(idx_cell))./(n_cntr(colIdx(fnd))-n_left(colIdx(fnd))).*(nvel(fnd)-n_left(colIdx(fnd))')';
        fnd=fright & inrange;
        idx_cell=sub2ind(size(msh(ct).Z),max(1,frw(colIdx(fnd))),colIdx(fnd)');
        vel_minz(fnd)=ZMmin(idx_cell)+(ZRmin(idx_cell)-ZMmin(idx_cell))./(n_right(colIdx(fnd))-n_cntr(colIdx(fnd))).*(nvel(fnd)-n_cntr(colIdx(fnd))')';
        f_incs=czvel<=maxz & czvel >= vel_minz & sigVel > sigmin;
        dZ=nan(size(nvel)); % Step in sigma
        fnd=f_incs & fleft & inrange; % Select all data in current section, whithin the section and on the left side of cells
        dZ(fnd)=dzed_left(colIdx(fnd))+(dzed_cnt(colIdx(fnd))-dzed_left(colIdx(fnd)))./(n_cntr(colIdx(fnd))-n_left(colIdx(fnd)))'.*(nvel(fnd)-n_left(colIdx(fnd))'); % Determine local vertical cell size (in sigma coordinates) for velocity estimates
        fnd=f_incs & fright & inrange; % Select all data in current section, whithin the section and on the right side 
        dZ(fnd)=dzed_cnt(colIdx(fnd))+(dzed_right(colIdx(fnd))-dzed_cnt(colIdx(fnd)))./(n_right(colIdx(fnd))-n_cntr(colIdx(fnd)))'.*(nvel(fnd)-n_cntr(colIdx(fnd))');% Determine local vertical cell size (in sigma coordinates) for velocity estimates 
        fnd=f_incs & inrange; % find velocity data belonging to current section, within the section
        rowIdx(fnd)=floor((maxz-czvel(fnd))./dZ(fnd))+1; % Calculate Index mapping to the right vertical cell       
    end
    
    if IsConventional
        rowIdx=repmat(rowIdx,[1 1 4]);
        colIdx=repmat(colIdx,[1 1 4]);
        f_incs=repmat(f_incs,[1 1 4]);
    end
    
    %% Velocity processing
    cMod=cat(4,vel(:,fcur,:),...                       % beam velocity (in conventional processing this is velocity in earth coordinates)
               bsxfun(@times,TM1(:,fcur,:),mod_ut(:,fcur,:,:)),... % model for u
               bsxfun(@times,TM2(:,fcur,:),mod_vt(:,fcur,:,:)),...;% model for v
               bsxfun(@times,TM3(:,fcur,:),mod_wt(:,fcur,:,:)) );   % model for w
    siz=[size(msh(ct).Z,1) size(msh(ct).Z,2)];
    datacol=cell([siz nanmax(tid(ct,:))]); % Initialize variable to collect velocity and tr. matrix data

    for ccr=nanmin(tid(ct,:)):nanmax(tid(ct,:)) % loop over all crossings
        if IsCumulative
            fnd=bsxfun(@and,tid(ct,fcur)<=ccr, f_incs & rowIdx<=size(msh(ct).Z,1) & rowIdx>0); % find data in current section, in all crossing up to current, and belonging to any cell
        else
            fnd=bsxfun(@and,tid(ct,fcur)==ccr, f_incs & rowIdx<=size(msh(ct).Z,1) & rowIdx>0); % find data in current section, in all crossing up to current, and belonging to any cell
        end
        if proxim>0 % if a proximity is given
            fnd=fnd & abs(svel)<=proxim; % only include data within proxim distance from the section
        end
        if ~any(fnd(:)), continue, end;
        col_dat=cell(npars+1,1);
        for cPars=1:npars+1
            select_par=false([size(fnd), npars+1]);
            select_par(:,:,:,cPars)=true;
            col_dat{cPars}=accumarray({rowIdx(fnd),colIdx(fnd)},cMod(bsxfun(@and,fnd,select_par)),siz,@(x) {x},{double.empty(0,1)});
        end
        datacol(:,:,ccr)=cellfun(@horzcat,col_dat{:},'UniformOutput',false);
    end
    [msh(ct).pars, msh(ct).rsq, msh(ct).S, msh(ct).Nvel]=cellfun(@(x) estvel(x,rmresid), datacol,'UniformOutput',false); % Estimate velocity with least squares estimates (see function estvel.m for procedure)
    msh(ct).pars=cell2mat(msh(ct).pars); % reshape velocity data
    msh(ct).rsq=cell2mat(msh(ct).rsq); % reshape mean squared error data
    msh(ct).Nvel=cell2mat(msh(ct).Nvel); % reshape Number of valid samples data
    msh(ct).S=cell2mat(msh(ct).S);

      
    %% Remove data with high std
    if nstd>0
        S=reshape(msh(ct).S(msh(ct).p.progfgood_tensdiag),[1 size(msh(ct).p.Z,2) n_rep_trans npars]);
        medS=nanmedian(S,2);
        fbad=any(bsxfun(@gt, S, medS*nstd^2),4);
        msh(ct).S(msh(ct).p.progfgood_tens((repmat(fbad(:),npars^2,1))))=nan;
        msh(ct).pars(msh(ct).p.progfgood_vec((repmat(fbad(:),npars,1))))=nan;
        msh(ct).Nvel(msh(ct).p.progfgood(fbad(:)))=nan;
    end
    
    %% Perform rotations
    if rotate
        if ~isequal(mod_ut,mod_vt)
            warning('procTrans:NotEqualUVMods','To perform a rotation the model for u and v should be identical, skipping rotations!')
        else
            npars_uv=npars_u;        
            % search for average velocity (i.e. model ==1) instead of assuming
            % it is in front!
            u_par=find(all(mod_ut(:,fcur,:,:)==1,2));
            v_par=find(all(mod_vt(:,fcur,:,:)==1,2));

            if numel(u_par)~=1 || numel(v_par)~=1
                warning('procTrans2:NoMeanInModel','Could not find mean velocity in u and v models, skipping rotations based on horizontal velocity')
            else
                % Rotation to cross-section averaged flow direction
                msh(ct).cs.dir=atan2(nanmean(nanmean(msh(ct).pars(:,:,:,npars_u+v_par),1),2),nanmean(nanmean(msh(ct).pars(:,:,:,u_par),1),2)); % Transect averaged flow direction
                cs=cos(msh(ct).cs.dir);
                ss=sin(msh(ct).cs.dir);
                csrot=zeros(1,1,nrepeat,npars,npars);
                idxrw=repmat([1:npars, 1:npars_uv, npars_uv+1:2*npars_uv],1,nrepeat);
                idxcol=repmat([1:npars, npars_uv+1:2*npars_uv, 1:npars_uv],1,nrepeat);
                idxrep=reshape(bsxfun(@times,ones(npars+2*npars_uv,1),1:nrepeat),1,[]);
                idxz=ones(size(idxrep));
                idxn=ones(size(idxrep));
                val=reshape(squeeze([repmat(cs,1,npars_uv*2) ones(1,npars_w,nrepeat) repmat(ss,1,npars_uv) repmat(-ss,1,npars_uv)]),1,[]);
                csrot(sub2ind(size(csrot),idxz,idxn,idxrep,idxrw,idxcol))=val;        
                st1=matmult(msh(ct).S,permute(csrot,[1 2 3 5 4]),[4 5]);   % Tensor rotation SR' (step 1)
                msh(ct).cs.S=matmult(csrot,st1,[4 5]); % Tensor rotation R(SR') (step 2)
                msh(ct).cs.pars=matmult(csrot,msh(ct).pars,[4 5]); % Vector rotation Rv

            %     % Vertically averaged flow direction
                msh(ct).da.dir=atan2(nanmean(msh(ct).pars(:,:,:,npars_u+v_par),1),nanmean(msh(ct).pars(:,:,:,u_par),1)); % depth averaged flow direction
                darot=zeros(1,siz(2),nrepeat,npars,npars);
                cs=cos(msh(ct).da.dir);
                ss=sin(msh(ct).da.dir);
                idxrw=repmat([1:npars, 1:npars_uv, npars_uv+1:2*npars_uv],1,nrepeat*siz(2));
                idxcol=repmat([1:npars, npars_uv+1:2*npars_uv, 1:npars_uv],1,nrepeat*siz(2));
                idxn=repmat(reshape(bsxfun(@times,ones((npars+2*npars_uv),1),1:siz(2)),1,[]),1,nrepeat);
                idxrep=reshape(bsxfun(@times,ones((npars+2*npars_uv)*siz(2),1),1:nrepeat),1,[]);
                idxz=ones(size(idxrep));     
                val=reshape(squeeze([repmat(cs,npars_uv*2,1,1); ones(npars_w,siz(2),nrepeat); repmat(ss,npars_uv,1,1); repmat(-ss,npars_uv,1,1)]),1,[]);
                darot(sub2ind(size(darot),idxz,idxn,idxrep,idxrw,idxcol))=val;
                st1=matmult(msh(ct).S,permute(darot,[1 2 3 5 4]),[4 5]);% Tensor rotation SR' (step 1)
                msh(ct).da.S=matmult(darot,st1,[4 5]);% Tensor rotation R(SR') (step 2)
                msh(ct).da.pars=matmult(darot,msh(ct).pars,[4 5]);% Vector rotation Rv
            end % if numel(u_par)~=1 || numel(v_par)~=1
        %     % Cross-section direction
            msh(ct).sec.dir=atan2(N(2),N(1));
            secrot=zeros(1,1,1,npars,npars);
            cs=cos(msh(ct).sec.dir);
            ss=sin(msh(ct).sec.dir);
            idxrw=[1:npars, 1:npars_uv, npars_uv+1:2*npars_uv];
            idxcol=[1:npars, npars_uv+1:2*npars_uv, 1:npars_uv];
            idxrep=ones(size(idxrw));
            idxz=ones(size(idxrep));
            idxn=ones(size(idxrep));
            val=[repmat(cs,1,npars_uv*2) ones(1,npars_w) repmat(ss,1,npars_uv) repmat(-ss,1,npars_uv)];
            secrot(sub2ind(size(secrot),idxz,idxn,idxrep,idxrw,idxcol))=val;
            st1=matmult(msh(ct).S,permute(secrot,[1 2 3 5 4]),[4 5]);   % Tensor rotation SR' (step 1)
            msh(ct).sec.S=matmult(secrot,st1,[4 5]);% Tensor rotation R(SR') (step 2)
            msh(ct).sec.pars=matmult(secrot,msh(ct).pars,[4 5]);% Vector rotation Rv
        end % if isequal(mod_ut,mod_vt)
    end % if rotate
end


function [pars, rsq, S, nin]=estvel(in,rmresid)
% estimates cartesian velocities given beam velocities and correspoding
% transformation matrix terms. IN is Nx4, N being the amount of beam
% velocity samples. Each row is composed of a beam velocity and the three
% terms from the transformation matrix

in(any(isnan(in),2),:)=[]; % Throw out bad data
rmvel=true; % initialize variable holding indices to outliers
n_pars=size(in,2)-1;
while any(rmvel) % iterate until no outliers are found
    if isempty(in) || rank(in(:,2:end))<n_pars % if no velocity is available or matrix is rank deficient
        pars=nan(n_pars,1); rsq=nan; S=nan(n_pars,n_pars); % return nans
        break; % Exit loop
    else
        [pars, ~, ~, S]=lscov(in(:,2:end),in(:,1)); % Perform least squares estimate of velocities
        res=abs(in(:,1)-in(:,2:end)*pars); % Determine residuals in beam velocity
        rsq=1-sum(res.^2)./sum((in(:,1)-mean(in(:,1))).^2);
    end
    if rmresid==0, break, end % if rmvel == 0, no outlier removal needed, exit loop
    rmvel=res>rmresid*median(res); % Detect outliers from residuals
    in(rmvel,:)=[]; % Remove outliers
end
nin=size(in,1); % compute amount of points used for final velocity computation
pars=shiftdim(pars,-3); % shift dimension to output 1x1x1x1x3 vector
rsq=shiftdim(rsq,-3); % shift dimension to output 1x1x1x1x3 vector
S=shiftdim(S,-3); % shift dimension to output 1x1x1x1x3x3 vector



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
