function msh=procTrans(adcp,tid,varargin)
% PROCTRANS processes repeat transect vessel mounted ADCP data
%
%   MSH=PROCTRANS(ADCP,TID) is processing the adcp-structure ADCP (as 
%       generated by readadcp or readDeployment). The input TID is an CSxE
%       matrix where CS is the number of cross-sections one whishes to
%       analyse and E is the number of ensembles. In each row a zero
%       indicates the ensemble is not part of the cross-section. All
%       ensembles with the same number in a row of TID are processed 
%       together. The output structure MSH is a CSx1 structure with the
%       following fields:
%       eta - average eta for each time step
%       maxeta - maximum eta for each time step
%       mineta - minimum eta for each time step
%       time - average time for each time step
%       Tvec - unit vector tangential to the cross-section in x-y plane
%       Nvec - unit vector orthogonal to the cross-section in x-y plane
%       Pm - center of the cross-section
%       N - RxC matrix with N coordinates (i.e. along Tvec) of cell-centers
%           with R the maximum number of rows and C the number of columns
%           in the mesh
%       X - X coordinates of cell-centers (RxC)
%       Y - RxC matrix with Y coordinates of cell-centers
%       Z - RxCxT array with Z coordinates of cell-centers for each 
%           time-step, with T the number of time-steps
%       Sig - RxCxT array with Sigma coordinates of cell-centers for each
%             time-step
%       pars - RxCxTxP array with estimated parameters. P is the number of
%              parameters
%       rsq - RxCxT array with the coefficient of determination of the 
%             parameter fit
%       Nvel - RxCxT array with the number of raw-beam velocities included
%             in the fit
%       S - RxCxTxPxP array with the covariance matrix (PxP) of the
%           estimated parameters.
%       p - structure to ease plotting of the results:
%           p.nbed - n coordinate of the bed (defined at cell sides and
%                    centers)
%           p.xbed - x coordinate of the bed
%           p.ybed - y coordinate of the bed
%           p.zbed - z coordinate of the bed
%           p.fgood - RxC logical matrix with 'true' for cells in the
%                     cross-section
%           p.progfgood - QxTx1 vector with the linear indices indicating 
%                         the cells in the cross-section for an RxCxT 
%                         matrix. Q is the number of cells within a 
%                         cross-section Allows to easily reference all 
%                         cells within the cross-section at a certain time.
%           p.progfgood_vec - QxTxP array with linear indices to reference
%                             a parameter for an RxCxTxP output array
%                             (pars)
%           p.progfgood_tens - QxTxPxP array with linear indices to
%                             reference a value for and RxCxTxPxP output
%                             array (S)
%           p.progfgood_tensdiag - QxTxP array with linear indices to
%                                  reference a diagonal element of a tensor
%                                  output RxCxTxPxP (S)
%           p.N - N coordinates of the cell vertices. Usefull when plotting
%                 the mesh with the PATCH command
%           p.X - X coordinates of the cell vertices (7xQ)
%           p.Y - Y coordinates of the cell vertices (7xQ)
%           p.Z - Z coordinates of the cell vertices (7xQxT)
%           p.Sig - Sigma coordinate of the cell vertices (7xQxT)
%       sec - Structure containing the fitting results rotated to match with
%            the direction of the cross-section (Nvec):
%            sec.dir - Rotation direction (mathematical angle, 1x1xT)
%            sec.pars - Rotated parameters (RxCxTxP)
%            sec.S - Rotated covariance matrix of parameters (RxCxTxP)
%       da - Structure containting fitting results rotated to the depth
%            averaged flow direction:
%            da.dir - Rotation direction (mathematical angle,  1xCxT)
%            da.pars - Rotated parameters (RxCxTxP)
%            da.S - Rotated covariance matrix of parameters (RxCxTxP)
%       cs - Structure containing fitting results rotated to the
%             cross-sectionally averaged flow direction
%             cs.dir - Rotation direction (mathematical angle, 1x1xT)
%             cs.pars - Rotated parameters (RxCxTxP)
%             cs.S - Rotated covariance matrix of parameters (RxCxTxP)
%       In the description above the following letters are used to
%       represent a certain dimension
%           CS - number of cross-section
%           E - number of ensembles
%           R - maximum number of rows in output mesh
%           C - number of columns in output mesh
%           T - number of time steps in output
%           P - number of fitted parameter (i.e. sum of the parameters for
%               x,y and z components of the velocity)
%           Q - number of cells in the cross-section
%
%    MSH=PROCTRANS(ADCP,TID,'ParamName', ParamValue) Allows to specify
%    several options for the processing:
%         'DepthTransducer'
%           scalar value - Sets the depth of the transducers in m. Default
%           is 0.3m.
%         'DeltaN'
%           scalar value - Sets the mesh size in n-direction. Default is 5m
%         'DeltaZ'
%           scalar value - Sets the mesh size in z-direction. Default is 1m
%         'Eta'
%           1xE row vector - Gives the value of eta (Water level) for each
%           ensmeble. Default is 0 (no water level changes)
%         'MinimumSigma'
%           scalar value between 0 and 1 - Sets the minium sigma to start
%           the mesh. Default is 0.06.
%         'CumulateCrossings'
%           true | {false} - Whether the time-steps given in TID for each
%           row should be cumulated
%         'ConventionalProcessing'
%           true | {false} - Whether to use the conventional processing.
%         'TopMeshLowestEta'
%           true | {false} - Whether to set the top of the mesh to the
%           lowest value of eta or to the top of the data thoughout the
%           dataset
%         'ConstantZetaMesh'
%           true | {false} - Whether to use a mesh constant in Z
%           coordinates instead of constant in Sigma coordinates
%         'RemoveOutliers'
%           positive scalar - Indicates how many times the residual in
%           beam-velocity should exceed the median of all residuals to 
%           discard a beam-velocity from the fit. Default is 0, indicating 
%           no outlier removal
%         'Proximity'
%           positive scalar - How close should the beam-velocity be in 
%           s-direction (i.e. along Nvec) to be included in the fit. 
%           Default is 0, indicating to include all data
%         'StdFiltering'
%           positive scalar - Indicates how many times the standard
%           deviation of any parameter in a cell should excess the median
%           of the standard deviations over all cells to be treated as bad.
%           Default is 6. A value of 0 deactivates this filtering.
%         'ShipReference'
%           {'bt'} | 'gps' | 'btgps' | 'gpsbt' - Indicates which method
%           should be used to compute the ship velocity:
%           'bt' use bottom-tracking
%           'gps' use GPS
%           'btgps' use bottom-tracking and GPS if bottom-tracking is
%               unavailable
%           'gpsbt' use GPS and bottom-tracking when GPS is unavailable
%         'Pusr'
%           two element row vector - Sets a user defined center of the
%           cross-section. Can be usefull when the measurements are not
%           uniformly collected along the cross-section
%         'RotatePars'
%           {true} | false - Whether to perform rotations
%         'Model'
%           scalar function handle - provides a function to specify a
%           velocity model. The function should accept as input nx5 matrix
%           where each column represents the n, z, sigma, s and t 
%           coordinates of a beam velocity. The function must return three
%           matrices representing the models for u, v and w respectively.
%           The matrices have n rows and as many columns as parameters in
%           the model. Default function fits a constant mean velocity for
%           u, v and w at each time step.
%         'Known'
%           scalar function handle - provides a function to specify the
%           know terms in the velocity model. The function should accept as
%           input nx5 matrix where each column represents the n, z, sigma, 
%           s and t coordinates of a beam velocity. The function must 
%           return three vectors representing the known terms for u, v and 
%           w respectively. 
%         'GetVelocity'
%           scalar function handle - defines how to retrieve velocity from
%           the estimated parameters. It accepts and nxP matrix as input.
%           It returns three matrices with the estimated u,v and w
%           velocity.

%       TODO: make fitting function for:
%       fourier fitting
%       log-profile


%    Copyright 2013,2014,2015 Bart Vermeulen
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




%% Handle input
assert(isADCPstruct(adcp)); % Check first input is an ADCP structure
assert(size(tid,2)==size(adcp.VEL,2)); % Ensure size of second dimension for the transect id's matches the number of ensembles in the data

nens=size(adcp.VEL,2);

P=inputParser;
P.FunctionName='procTrans';
P.addParameter('DepthTransducer',0.3,@(x) isscalar(x) && isnumeric(x) && x>0) % depth of transducer
P.addParameter('DeltaN',5,@(x) isscalar(x) && isnumeric(x) && x>0) % N resolution of mesh
P.addParameter('DeltaZ',1,@(x) isscalar(x) && isnumeric(x) && x>0) % target Z resolution of mesh
P.addParameter('Eta',0,@(x) isvector(x) && isnumeric(x) && (numel(x)==1 || numel(x)==nens)) % eta (water level)
P.addParameter('MinimumSigma',0.06,@(x) isscalar(x) && isnumeric(x) && x>=0 && x <=1) % Minimum sigma for meshing (to account for side lobes)
P.addParameter('CumulateCrossings',false,@(x) islogical(x) && isscalar(x)) % Set whether to lump consecutive crossings (only works if crossings are indicated in tid)
P.addParameter('ConventionalProcessing',false,@(x) islogical(x) && isscalar(x)) % Set whether to use conventional processing
P.addParameter('TopMeshLowestEta',false,@(x) islogical(x) && isscalar(x)) % Set whether top of mesh should be set at lowest water level (or top of data if this is lower than lowest water level)
P.addParameter('ConstantZetaMesh',false,@(x) islogical(x) && isscalar(x)) % Set whether to keep sigma constant. If false z is kept constant
P.addParameter('RemoveOutliers',0,@(x) isscalar(x) && isnumeric(x) && x>=0); % Remove outliers when their residuals exceed this amount of times the median residual in velocity inversion
P.addParameter('Proximity',0,@(x) isscalar(x) && isnumeric(x) && x>=0); % Set the maximum distance from the cross-section for data to be included in the calculation
P.addParameter('StdFiltering',6,@(x) isscalar(x) && isnumeric(x) && x>=0); % Remove velocity with a standard deviation of the estimate exceeding this amount of times the median standard deviation
P.addParameter('ShipReference','bt',@(x) ischar(x) && any(strcmpi(x,{'bt','gps','gpsbt','btgps'}))); % Select which boat velocity calculation to use
P.addParameter('Pusr',[], @(x) isnumeric(x)); % User set the central position of a transect (can be usefull when data density is not uniform on transect
P.addParameter('RotatePars',true,@(x) isscalar(x) && islogical(x));
P.addParameter('EnableDebugging',false,@(x) isscalar(x) && islogical(x));
P.addParameter('Model',@default_model,@(x) isscalar(x) && isa(x,'function_handle'));
P.addParameter('Known',@default_known,@(x) isscalar(x) && isa(x,'function_handle'));
P.addParameter('GetVelocity',@default_getvelocity,@(x) isscalar(x) && isa(x,'function_handle'));

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

water_level=ConstantWaterLevel(eta);
eta=water_level.get_water_level(time);

shref=P.Results.ShipReference;
Pusr=P.Results.Pusr;
IsCumulative=P.Results.CumulateCrossings;
IsConventional=P.Results.ConventionalProcessing;
rotate=P.Results.RotatePars;
if rotate && ~(any(strcmp(P.UsingDefaults,'Model')) && any(strcmp(P.UsingDefaults,'Known'))) && any(strcmp(P.UsingDefaults,'GetVelocity'))
    rotate=false;
    warning('When defining a custom velocity model, rotations of velocity are only possible when also a velocity retrieval function is provided')
end
fdebug=P.Results.EnableDebugging;

f_model=P.Results.Model;
f_known=P.Results.Known;
f_velocity=P.Results.GetVelocity;

%% Get ADCP positioning data
[xadcp,yadcp]=utmADCP(adcp); % Get adcp position
misal=getExtMisalign(adcp); % Get beam 3 misalignment

% depth location
[xb,yb,zb]=depthADCP(adcp,'Beam3Misalign',misal); % Get offsets from ADCP to the bed
zb=bsxfun(@plus,zb-depthtransd,eta); % transform vertical offset to bed level with respect to eta=0
xb=bsxfun(@plus,xadcp',xb); % transform offsets to UTM coordinate system for horizontal
yb=bsxfun(@plus,yadcp',yb); % transform offsets to UTM coordinate system for horizontal
if IsConventional
    xb=repmat(nanmean(xb,3),[1 1 4]); % determine beam average for conventional velocity processing
    yb=repmat(nanmean(yb,3),[1 1 4]); % determine beam average for conventional velocity processing
    zb=repmat(nanmean(zb,3),[1 1 4]); % determine beam average for conventional velocity processing
end


% velocity location
[xvel,yvel,zvel]=mapADCP(adcp,'IsUpward', false, 'Beam3Misalign',misal); % get offsets from adcp to velocity locations
zvel=bsxfun(@plus,zvel-depthtransd,eta); % transform vertical offset to z position of velocity with respect to eta=0
xvel=bsxfun(@plus,xadcp',xvel); % transform offsets to vel data to UTM coordinate system for horizontal
yvel=bsxfun(@plus,yadcp',yvel); % transform offsets to vel data to UTM coordinate system for horizontal
if IsConventional
    xvel=repmat(nanmean(xvel,3),[1 1 4]); % determine beam average for conventional velocity processing
    yvel=repmat(nanmean(yvel,3),[1 1 4]); % determine beam average for conventional velocity processing
    zvel=repmat(nanmean(zvel,3),[1 1 4]); % determine beam average for conventional velocity processing
end
tvel=repmat(time,[size(xvel,1),1,size(xvel,3)]);

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
        [bvel,~]=corADCP(adcp2,'b','UseExtHeading',true,'Beam3Misalign',misal,'ForceOrigin','e'); % transform to earth velocity;
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
mesh(nsec)=Mesh;
bathy(nsec)=BathymetryScatteredPoints;
[bathy(1:nsec).water_level]=deal(water_level);

% %% create time index (mapping to right time cell)
% eta = eta .* ones(1,nens); % make sure eta is row vector with same number of elements as number of ensembles

% Start processing
for ct=1:size(tid,1) % For all sections  
    fcur=tid(ct,:)>0; % mask, indicating data belonging to current section
    ctvel=tvel(:,fcur,:);
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
        mesh(ct).origin = mean(Pusr(:,:,ct),2);
    else
        mesh(ct).origin=nanmean(P,2); % Average ADCP position vector
    end
    mesh(ct).direction=princdir(P); % Section tangential vector, determined as largest eigenvector of P
    
    % Project positions (inner product with unit vectors)
    posb=mesh(ct).xy2sn([xb(:,fcur,:);yb(:,fcur,:)]);
    posvel=mesh(ct).xy2sn([shiftdim(xvel(:,fcur,:),-1);shiftdim(yvel(:,fcur,:),-1)]);
    
    % Store important vectors
    msh(ct).Tvec=mesh(ct).direction; % tangential vector
    msh(ct).Nvec=mesh(ct).direction_orthogonal; % orthogonal vector
    msh(ct).Pm=mesh(ct).origin; % mean position vector
    
    %% Calculate sigma for velocity locations
    posb=[posb; zb(:,fcur,:)];   
    
    % Finite depth and velocity masks
    fgoodb=repmat(all(isfinite(posb),1),3,1); % Mask for finite bed positions
    fgoodvel=all(isfinite(posvel),1); % Masf for finite velocity positions

    % Initialize variables
    zbvel=nan(size(fgoodvel)); % Bed elevation at velocity locations
    
    bathy(ct).pos=reshape(posb(fgoodb),3,[]);
    zbvel(fgoodvel)=bathy(ct).get_bed_elev(reshape(posvel(repmat(fgoodvel,2,1)),2,[]));
%     posvel=[posvel;zbvel];
    
    % Compute sigma
    czvel=shiftdim(zvel(:,fcur,:),-1);
    posvel=[posvel; czvel];
    if numel(eta)>1
        ceta=eta(fcur); 
    else
        ceta=eta; 
    end
    ceta=shiftdim(ceta,-1);
    sigVel=nan(size(czvel));
    f_incs=bsxfun(@gt,czvel,zbvel);
    dVel=bsxfun(@minus,ceta,zbvel);
    sigVel(f_incs)=(czvel(f_incs)-zbvel(f_incs))./dVel(f_incs); % Compute sigma for velocity locations
    f_incs=sigVel<=1 & sigVel>=sigmin; % Create mask for velocity locations in the cross-section
    
    
    %% Meshing
    % Compute minimum and maximum n
    nmin=nanmin(posvel(f_incs & [false; true; false])); % determine minimum n
%     Pn=mesh(ct).origin+mesh(ct).direction*nmin; % vector pointing to minimum n
%     posvel=posvel-[0; nmin; 0]; % remove minimum n from velocity n coordinates (now they range between 0 and maximum n)
%     posb=posb-[0; nmin; 0]; % remove minimum n from depth n coordinates (now they range between 0 and maximum n)
    nmax=nanmax(posvel(f_incs & [false; true; false])); % compute maximum n
    num_verticals=ceil((nmax-nmin)./veldn); % determine number of vertical in section

    delta_n=(nmax-nmin)/num_verticals;
    n_cell_left=nmin+((1:num_verticals)-1)*delta_n;
    n_cell_center=n_cell_left+delta_n/2;
    n_cell_right=n_cell_left+delta_n;
%     xy_cell_left=mesh(ct).origin+n_cell_left.*mesh(ct).direction;
%     xy_cell_right=mesh(ct).origin+n_cell_left.*mesh(ct).direction;
%     xy_cell_center=mesh(ct).origin+n_cell_left.*mesh(ct).direction;
    zb_cell_left=bathy(ct).get_bed_elev(n_cell_left.*[0; 1]);
    zb_cell_right=bathy(ct).get_bed_elev(n_cell_right.*[0; 1]);
    zb_cell_center=bathy(ct).get_bed_elev(n_cell_center.*[0; 1]);
    
    

    % Compute maximum z verticals
    maxz=nanmax(czvel(:));
    minz_cell_left=zb_cell_left+(water_level.get_water_level()-zb_cell_left)*sigmin;
    minz_cell_center=zb_cell_center+(water_level.get_water_level()-zb_cell_center)*sigmin;
    minz_cell_right=zb_cell_right+(water_level.get_water_level()-zb_cell_right)*sigmin;

    num_cells_per_vert= ceil((maxz-minz_cell_center)./veldz);
    num_cells=sum(num_cells_per_vert);
    mesh(ct).cells(num_cells)=MeshCell;
    dz_center=(maxz-minz_cell_center)./num_cells_per_vert;
    dz_left=(maxz-minz_cell_left)./num_cells_per_vert;
    dz_right=(maxz-minz_cell_right)./num_cells_per_vert;
    max_num_cells=max(num_cells_per_vert);
    mask_mesh=(1:max_num_cells)'<=num_cells_per_vert';
    [col_number,row_number]=meshgrid(1:num_verticals,1:max_num_cells);
    col_idx=col_number(mask_mesh);
    row_idx=row_number(mask_mesh);
    nz_coords=[permute(cat(3, n_cell_left(col_idx), n_cell_center(col_idx), n_cell_right(col_idx), n_cell_right(col_idx), n_cell_center(col_idx), n_cell_left(col_idx), n_cell_left(col_idx)), [1 3 2]);...
            permute(cat(3,maxz-dz_left(col_idx).*row_idx, maxz-dz_center(col_idx).*row_idx, maxz-dz_right(col_idx).*row_idx, maxz-dz_right(col_idx).*(row_idx-1), maxz-dz_center(col_idx).*(row_idx-1), maxz-dz_left(col_idx).*(row_idx-1), maxz-dz_left(col_idx).*row_idx),[2 3 1])];
    mesh(ct).bed_position=[n_cell_left(1) n_cell_center n_cell_right(end); zb_cell_left(1) zb_cell_center' zb_cell_right(end)];
    sig_coords=(nz_coords(2,:,:)-permute(cat(3,zb_cell_left(col_idx), zb_cell_center(col_idx), zb_cell_right(col_idx), zb_cell_right(col_idx), zb_cell_center(col_idx), zb_cell_left(col_idx), zb_cell_left(col_idx)),[2 3 1]))./...
        (water_level.get_water_level() -  permute(cat(3,zb_cell_left(col_idx), zb_cell_center(col_idx), zb_cell_right(col_idx), zb_cell_right(col_idx), zb_cell_center(col_idx), zb_cell_left(col_idx), zb_cell_left(col_idx)),[2 3 1]));
    
    mesh(ct).coordinates=[nz_coords(1,:,:); sig_coords];
    
    %% Velocity processing
    siz=size(mask_mesh); % size of output mesh
    cvel=vel(:,fcur,:);
    cTM1=TM1(:,fcur,:);
    cTM2=TM2(:,fcur,:);
    cTM3=TM3(:,fcur,:);
    datacol=cell([siz nanmax(tid(ct,:))]); % Initialize variable to collect velocity and tr. matrix data

    for ccr=nanmin(tid(ct,tid(ct,:)>0)):nanmax(tid(ct,:)) % loop over all crossings
        fnd=bsxfun(@and,tid(ct,fcur)==ccr, f_incs & rowIdx<=size(msh(ct).Z,1) & rowIdx>0); % find data in current section, in current crossing, and belonging to any cell
        if proxim>0 % if a proximity is given
            fnd=fnd & abs(svel)<=proxim; % only include data within proxim distance from the section
        end
        if ~any(fnd(:)), continue, end;
            

        % compute model coefficients from dz, dn, ds and dt
        % Compute dz, dn, ds and dt and collect them for each mesh cell
        coldz=cellfun(@minus,...
            accumarray({rowIdx(fnd), colIdx(fnd)},czvel(fnd),siz,@(x) {x}, {double.empty(0,1)}),...
            num2cell(msh(ct).Z(:,:,ccr)),...
            'UniformOutput',false);
        colSig=cellfun(@minus,...
            accumarray({rowIdx(fnd), colIdx(fnd)},sigVel(fnd),siz,@(x) {x}, {double.empty(0,1)}),...
            num2cell(msh(ct).Z(:,:,ccr)),...
            'UniformOutput',false);
        coldn=cellfun(@minus,...
            accumarray({rowIdx(fnd), colIdx(fnd)},nvel(fnd),siz,@(x) {x}, {double.empty(0,1)}),...
            num2cell(msh(ct).N(:,:)),...
            'UniformOutput',false);
        colds=accumarray({rowIdx(fnd), colIdx(fnd)},svel(fnd),siz,@(x) {x}, {double.empty(0,1)});
        coldt=cellfun(@(x,y) (x-y)*24*3600,...
        accumarray({rowIdx(fnd), colIdx(fnd)},ctvel(fnd),siz,@(x) {x}, {double.empty(0,1)}),...
        repmat({msh(ct).time(ccr)},siz),...
        'UniformOutput',false);            
        [colmu,colmv,colmw]=cellfun(@(n,z,Sig,s,t) f_model([n,z,Sig,s,t]),coldn, coldz, colSig, colds,coldt,'UniformOutput',false);
        [colku,colkv,colkw]=cellfun(@(n,z,Sig,s,t) f_known([n,z,Sig,s,t]),coldn, coldz, colSig, colds,coldt,'UniformOutput',false);
        npars_u=size(colmu{1},2);
        npars_v=size(colmv{1},2);
        npars_w=size(colmv{1},2);
        npars=npars_u+npars_v+npars_w;
        if rotate && ~all(all(cellfun(@isequal,colmu,colmv)))
            rotate = false;
            warning('procTrans:NotEqualUVMods','To perform a rotation the model for u and v should be identical, skipping rotations!')        
        end
        % collect tilt transformation matrix terms
        colTM1=accumarray({rowIdx(fnd),colIdx(fnd)},cTM1(fnd),siz,@(x) {x},{double.empty(0,1)});
        colTM2=accumarray({rowIdx(fnd),colIdx(fnd)},cTM2(fnd),siz,@(x) {x},{double.empty(0,1)});
        colTM3=accumarray({rowIdx(fnd),colIdx(fnd)},cTM3(fnd),siz,@(x) {x},{double.empty(0,1)});
        % collect beam velocity
        colbvel=accumarray({rowIdx(fnd),colIdx(fnd)},cvel(fnd),siz,@(x) {x},{double.empty(0,1)});
        % compute all terms for velocity inversion
        datacol(:,:,ccr)=cellfun(@(b,t1,t2,t3,mu,mv,mw,ku,kv,kw) ...
                                [b-(bsxfun(@times,ku,t1)+bsxfun(@times,kv,t2)+bsxfun(@times,kw,t3)) ...
                                 bsxfun(@times,mu,t1) ...
                                 bsxfun(@times,mv,t2) ...
                                 bsxfun(@times,mw,t3)],...
                          colbvel, colTM1, colTM2, colTM3, colmu, colmv, colmw, colku, colkv, colkw,...
                          'UniformOutput',false);
    end
    [msh(ct).pars, msh(ct).rsq, msh(ct).S, msh(ct).Nvel]=cellfun(@(x) estvel(x,rmresid), datacol,'UniformOutput',false); % Estimate velocity with least squares estimates (see function estvel.m for procedure)
    msh(ct).pars=cell2mat(msh(ct).pars); % reshape velocity data
    msh(ct).rsq=cell2mat(msh(ct).rsq); % reshape mean squared error data
    msh(ct).Nvel=cell2mat(msh(ct).Nvel); % reshape Number of valid samples data
    msh(ct).S=cell2mat(msh(ct).S); % reshape covariance matrices

    %% Generate fgood for vector and tensor data
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
   
    %% Remove data with high std
    if nstd>0
        S=reshape(msh(ct).S(msh(ct).p.progfgood_tensdiag),[1 size(msh(ct).p.Z,2) n_rep_trans npars]);
        medS=nanmedian(S,2);
        fbad=any(bsxfun(@gt, S, medS*nstd^2),4);
        msh(ct).S(msh(ct).p.progfgood_tens((repmat(fbad(:),[npars^2,1]))))=nan;
        msh(ct).pars(msh(ct).p.progfgood_vec((repmat(fbad(:),[npars,1]))))=nan;
        msh(ct).Nvel(msh(ct).p.progfgood(fbad(:)))=nan;
    end
    
    
    %% Perform rotations
    
    if npars<1
        warning('No parameters estimated, not performing rotations')
        rotate=false;
    end
        
    if rotate
        npars_uv=npars_u;
        [u, v]=f_velocity(reshape(msh(ct).pars,[],npars));
        sizvel=size(msh(ct).pars); sizvel(end)=1;
        u=reshape(u,sizvel);
        v=reshape(v,sizvel);
        msh(ct).cs.dir=atan2(nanmean(nanmean(v,1),2),nanmean(nanmean(u,1),2)); % Transect averaged flow direction
        cs=cos(msh(ct).cs.dir);
        ss=sin(msh(ct).cs.dir);
        csrot=zeros(1,1,nrepeat,npars,npars);
        idxrw=repmat([1:npars, 1:npars_uv, npars_uv+1:2*npars_uv],[1,nrepeat]);
        idxcol=repmat([1:npars, npars_uv+1:2*npars_uv, 1:npars_uv],[1,nrepeat]);
        idxrep=reshape(bsxfun(@times,ones(npars+2*npars_uv,1),1:nrepeat),1,[]);
        idxz=ones(size(idxrep));
        idxn=ones(size(idxrep));
        val=reshape(squeeze([repmat(cs,[1,npars_uv*2]) ones(1,npars_w,nrepeat) repmat(ss,[1,npars_uv]) repmat(-ss,[1,npars_uv])]),1,[]);
        csrot(sub2ind(size(csrot),idxz,idxn,idxrep,idxrw,idxcol))=val;        
        st1=matmult(msh(ct).S,permute(csrot,[1 2 3 5 4]),[4 5]);   % Tensor rotation SR' (step 1)
        msh(ct).cs.S=matmult(csrot,st1,[4 5]); % Tensor rotation R(SR') (step 2)
        msh(ct).cs.pars=matmult(csrot,msh(ct).pars,[4 5]); % Vector rotation Rv

        % Vertically averaged flow direction
        msh(ct).da.dir=atan2(nanmean(v,1),nanmean(u,1)); % depth averaged flow direction
        darot=zeros(1,siz(2),nrepeat,npars,npars);
        cs=cos(msh(ct).da.dir);
        ss=sin(msh(ct).da.dir);
        idxrw=repmat([1:npars, 1:npars_uv, npars_uv+1:2*npars_uv],[1,nrepeat*siz(2)]);
        idxcol=repmat([1:npars, npars_uv+1:2*npars_uv, 1:npars_uv],[1,nrepeat*siz(2)]);
        idxn=repmat(reshape(bsxfun(@times,ones((npars+2*npars_uv),1),1:siz(2)),1,[]),[1,nrepeat]);
        idxrep=reshape(bsxfun(@times,ones((npars+2*npars_uv)*siz(2),1),1:nrepeat),1,[]);
        idxz=ones(size(idxrep));     
        val=reshape(squeeze([repmat(cs,[npars_uv*2,1,1]); ones(npars_w,siz(2),nrepeat); repmat(ss,[npars_uv,1,1]); repmat(-ss,[npars_uv,1,1])]),1,[]);
        darot(sub2ind(size(darot),idxz,idxn,idxrep,idxrw,idxcol))=val;
        st1=matmult(msh(ct).S,permute(darot,[1 2 3 5 4]),[4 5]);% Tensor rotation SR' (step 1)
        msh(ct).da.S=matmult(darot,st1,[4 5]);% Tensor rotation R(SR') (step 2)
        msh(ct).da.pars=matmult(darot,msh(ct).pars,[4 5]);% Vector rotation Rv

        % Cross-section direction
        msh(ct).sec.dir=atan2(N(2,ct),N(1,ct));
        secrot=zeros(1,1,1,npars,npars);
        cs=cos(msh(ct).sec.dir);
        ss=sin(msh(ct).sec.dir);
        idxrw=[1:npars, 1:npars_uv, npars_uv+1:2*npars_uv];
        idxcol=[1:npars, npars_uv+1:2*npars_uv, 1:npars_uv];
        idxrep=ones(size(idxrw));
        idxz=ones(size(idxrep));
        idxn=ones(size(idxrep));
        val=[repmat(cs,[1,npars_uv*2]) ones(1,npars_w) repmat(ss,[1,npars_uv]) repmat(-ss,[1,npars_uv])];
        secrot(sub2ind(size(secrot),idxz,idxn,idxrep,idxrw,idxcol))=val;
        st1=matmult(msh(ct).S,permute(secrot,[1 2 3 5 4]),[4 5]);   % Tensor rotation SR' (step 1)
        msh(ct).sec.S=matmult(secrot,st1,[4 5]);% Tensor rotation R(SR') (step 2)
        msh(ct).sec.pars=matmult(secrot,msh(ct).pars,[4 5]);% Vector rotation Rv
    end % if rotate
end % for ct


end % function

function [mu, mv, mw]=default_model(in)
    [mu,mv,mw]=deal(ones(size(in,1),1));
end
function [ku, kv, kw]=default_known(in)
    [ku,kv,kw]=deal(zeros(size(in,1),1));
end
function [u, v, w]=default_getvelocity(in)
    u=in(:,1);
    v=in(:,2);
    w=in(:,3);
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
    if isempty(in) || n_pars < 1 || rank(in(:,2:end))<n_pars % if no velocity is available or matrix is rank deficient
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

end

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
end