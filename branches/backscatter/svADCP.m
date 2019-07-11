function SV = svADCP(inadcp,varargin)
% svADCP compute signal backscatter strenght (dB)
%        SV = svADCP(ADCP) compute signal backscatter strength from ADCP
%        profile data with a working version of the sonar equation (Deines,
%        1999):
% 
%        Sv = C + 10log10((Tx + 273.16)R2) - LDBM - PDBW + 2alphaR + Kc(E - Er)
% 
%        where
%
%        Sv is the backscattering strenght in dB referenced to 1/(4Pi m)
%        LDBM is 10log10(transmit pulse length, m)
%        PDBW is 10log10(transmit power, W)
%        Tx is temperature of the transducer (C)
%        R is range along the beam (slant range) to the scatterers (m)
%        alpha is absorption coefficient of water (dB/m)
%        C a constant +/- 3dB
%
%        SV = svADCP(ADCP,EnsRange) compute backscatter strenght for
%        selected ensemble range, EnsRange is a two column vector with
%        starting and ending indexes
%        SV = svADCP(ADCP,'isHADCP',true) compute backscatter strenght for
%        HADCP system
%
% See also: 
%


% TODO: add parameter for reference level Er
% TODO: add parameter for Kc
% TODO: add parameter for t_offset
% TODO: transducer diameter computation

%% Parsing input
P = inputParser;                                                           % Create an input parser
P.FunctionName = 'svADCP';                                                 % Check for correct spelling in function
P.addRequired('inadcp',@isstruct);                                         % Check adcp structure
P.addOptional('EnsRange',[0 0],@(x) isnumeric(x) && numel(x)==2);          % Ensemble range to process
P.addParameter('IsHadcp',false,@(x) islogical(x) && isscalar(x));          % isHADCP
% P.addParameter('alpha',0.1389,@(x) isscalar(x) && isnumeric(x));           % water attenuation

P.parse(inadcp,varargin{:});                                               % Parse input

EnsRange = P.Results.EnsRange;
IsHadcp = P.Results.IsHadcp;
% alpha=P.Results.alpha;
Kc = 0.45; % defaul from Deines                                            % RSSI scale factor dB per count (calibration !) Note that Aquavision uses 0.43
Er = 40;                                                                   % RSSI in a bucket, see PT3 command
t_offset = -0.35 ;                                                        % Temp Sens Offset (PS0 results)


clear P

%% Initialize parameters
nens = numel(inadcp.ensnum);                                               %Find number of ensembles
nbins = max(double(inadcp.nbins));                                              %Get the number of bins
if isequal(EnsRange,[0 0])
    EnsRange = [1 nens];
end
EnsIndex = EnsRange(1):EnsRange(end);
nens = range(EnsIndex) + 1;

%% Step 1: obtain ADCP characteristics measured at RDI's factory
bangle=get_beam_angle(inadcp);
pt=acoustics.PistonTransducer;
pt.water.salinity=35000;
% constants from rdi manuals
%TODO: Find correct transducer radii
current_fact=11451/1e6;
switch inadcp.sysconf(1,1:3) 
    case '000'
        pt.frequency=76.8e3;
        volt_fact=2092719/1e6;
        current_fact=43838/1e6;
        C=-159.1;
    case '100'
        pt.frequency=153.6e3;
        volt_fact=592157/1e6;
        C=-153.0;
    case '010'
        pt.frequency=307.2e3;
        volt_fact=592157/1e6;
        C=-143;
    case '110'
        pt.frequency=614.4e3;
        pt.radius=0.0505;
        volt_fact=380667/1e6;
        C=-139.3;
    case '001'
        pt.frequency=1228.8e3;
        pt.radius=0.027; % CHECK THIS
        volt_fact=253765/1e6;
        C=-129.1;
    case '101'
        pt.frequency=2457.6e3;
        pt.radius=0.0125; % CHECK THIS
        volt_fact=253765/1e6;
        C=NaN;
end
nbeams = 4;                                                                % number of beams


%% Step 3: obtain selected inadcp parameters
B = double(inadcp.blnk(inadcp.FileNumber(EnsIndex)))/100;                  % blank after transmit (m)
L = double(inadcp.lngthtranspulse(inadcp.FileNumber(EnsIndex)))/100;       % transmit pulse length (m)
LDBM = 10*log10(L);                                                         
D = double(inadcp.binsize(inadcp.FileNumber(EnsIndex)))/100;               % depth cell length (m)
curr = current_fact*double(inadcp.ADC(EnsIndex,1));                            % current (A)
volt = volt_fact*double(inadcp.ADC(EnsIndex,2));                            % voltage (V)
% Are these dependent on instrument ? Then, check for the HADCP !!!
DC_COEF = 9.82697464e1;                                                    % Temperature coefficients
FIRST_COEF = -5.86074151382e-3;
SECOND_COEF = 1.60433886495e-7;
THIRD_COEF = -2.32924716883e-12;
t_cnts = double(inadcp.ADC(EnsIndex,6))*256;                               % Temperature Counts (ADC value)
Tx = t_offset + ((THIRD_COEF.*t_cnts + SECOND_COEF).*t_cnts + FIRST_COEF).*t_cnts + DC_COEF; % real-time temperature of the transducer (C)

%% Step 4: prescribe relevant external variables
% sound absorption coefficient for each depth cell (dB/m)
% alpha = afromCTD(CTD);
% c = double(inadcp.speedsound(EnsIndex)) ;                                % speed of sound for each ensemble (m/s)
% creal = cfromCTD(CTD);

%% Step 5: determine transmit power
power = volt.*curr;                                                        % Transmit power in Watts
PDBW = 10*log10(power);

%% Step 6: evaluate variables in sonar equation for each depth cell
R = (repmat(B + (D+L)/2,nbins,1) + repmat((1:nbins)'-1,1,nens).*repmat(D,nbins,1) + repmat(D,nbins,1)/4)./repmat(cosd(bangle),nbins,1);                     %*creal/c  if correction for speed of sound is needed
                                                                           % creal is the average sound speed from transducer to the range cell.
% Rcritical = pi*Ray_dist/4;                                                 % R must be greater than Rcritical
% belowRcrit = R < Rcritical;                                                % bins below Rcritical
alpha=pt.attenuation;
alpha_n = 2*alpha*D./cosd(bangle);                                           % absorption for each range cell
two_alpha_R = 2*alpha*repmat(B,nbins,1)./repmat(cosd(bangle),nbins,1) + repmat((1:nbins)',1,nens).*repmat(alpha_n,nbins,1);                   % compute 2alphaR

%% Step 7: calculate profiles of backscatter coefficient
Tx = repmat(Tx(:)',[nbins 1 nbeams]);                                      % replicate matrices for computation
PDBW = repmat(PDBW(:)',[nbins 1 nbeams]);
LDBM = repmat(LDBM(:)',[nbins 1 nbeams]);
R = repmat(R,[1 1 nbeams]);
two_alpha_R = repmat(two_alpha_R,[1 1 nbeams]);
ECHO = double(inadcp.ECHO(:,EnsIndex,:));                                  % echo intensity (counts)

SV = C + 10*log10((Tx+273.16).*R.^2.*pt.near_field_correction(R).^2) - LDBM - PDBW + two_alpha_R + Kc.*(ECHO-Er);

end