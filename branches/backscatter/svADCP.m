function SV = svADCP(inadcp,varargin)
% svADCP compute signal backscatter strenght (dB)
%        SV = svADCP(ADCP) compute signal backscatter strenght from ADCP
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
% Last edit: 14/07/2009 not yet added Near Field correction


%% Parsing input
P = inputParser;                                                           % Create an input parser
P.FunctionName = 'svADCP';                                                 % Check for correct spelling in function
P.addRequired('inadcp',@isstruct);                                         % Check adcp structure
P.addOptional('EnsRange',[0 0],@(x) isnumeric(x) && numel(x)==2);          % Ensemble range to process
P.addParamValue('IsHadcp',false,@(x) islogical(x) && isscalar(x));         % isHADCP
P.addParamValue('alpha',0.1389,@(x) isscalar(x) && isnumeric(x));         % isHADCP

P.parse(inadcp,varargin{:});                                               % Parse input

EnsRange = P.Results.EnsRange;
IsHadcp = P.Results.IsHadcp;
alpha=P.Results.alpha;
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
if IsHadcp
    bangle = double(inadcp.HADCPbeamangle(inadcp.FileNumber(EnsIndex)));                                % determine beam angle for HADCP
    C = -139.3;                                                            % dB (everything we cannot measure)
    %PDBW = 9;                                                              % Expected Transmit Power (Watts)
    %Ray_dist = 1.96;                                                       % Rayleigh distance (m)
else
    bangle=bin2dec(inadcp.sysconf(:,9:10))';
    bangle(bangle==0)=15;
    bangle(bangle==2)=20;
    bangle(bangle==3)=30;
    bangle=bangle(inadcp.FileNumber(EnsIndex));
    C = -129.1;                                                            % dB (everything we cannot measure)
    %PDBW = 4.8;                                                            % Expected Transmit Power (Watts)
    %Ray_dist = 1.67;                                                       % Rayleigh distance (m)
end
if bangle == 0;
    warning('filterADCP2:UnknownBangle','Beam angle unclear, assuming 20 degrees')
    bangle = 20;
end
% switch inadcp.sysconf(1:3)
%     case '000'
%         NF=2*3.15; % 75 KHz
%     case '100'
%         NF=2*2.19; % 150 KHz
%     case '010'
%         NF=2*2.84; % 300 KHz
%     case '110'
%         NF=2*3.28; % 600 KHz
%     case '001'
%         NF=2*1.88; % 1200 KHz
% end
nbeams = 4;                                                                % number of beams
Kc = 0.45;                                                                 % RSSI scale factor dB per count (calibration !)

%% Step 2: calibrate reference level for echo intensity (PT3 results)
Er = 40;                                                                   % RSSI in a bucket

%% Step 3: obtain selected inadcp parameters
B = double(inadcp.blnk(inadcp.FileNumber(EnsIndex)))/100;                                               % blank after transmit (m)
L = double(inadcp.lngthtranspulse(inadcp.FileNumber(EnsIndex)))/100;                                    % transmit pulse length (m)
LDBM = 10*log10(L);
D = double(inadcp.binsize(inadcp.FileNumber(EnsIndex)))/100;                                            % depth cell length (m)
curr = 0.011451*double(inadcp.ADC(EnsIndex,1));                            % current (A)
volt = 0.253765*double(inadcp.ADC(EnsIndex,2));                            % voltage (V)
% Are these dependent on instrument ? Then, check for the HADCP !!!
DC_COEF = 9.82697464e1;                                                    % Temperature coefficients
FIRST_COEF = -5.86074151382e-3;
SECOND_COEF = 1.60433886495e-7;
THIRD_COEF = -2.32924716883e-12;
t_cnts = double(inadcp.ADC(EnsIndex,6))*256;                               % Temperature Counts (ADC value)
t_offset = -0.35 ;                                                         % Temp Sens Offset (PS0 results)
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
alpha_n = 2*alpha*D./cosd(bangle);                                           % absorption for each range cell
two_alpha_R = 2*alpha*repmat(B,nbins,1)./repmat(cosd(bangle),nbins,1) + repmat((1:nbins)',1,nens).*repmat(alpha_n,nbins,1);                   % compute 2alphaR

%% Step 7: calculate profiles of backscatter coefficient
Tx = repmat(Tx(:)',[nbins 1 nbeams]);                                      % replicate matrices for computation
PDBW = repmat(PDBW(:)',[nbins 1 nbeams]);
LDBM = repmat(LDBM(:)',[nbins 1 nbeams]);
R = repmat(R,[1 1 nbeams]);
two_alpha_R = repmat(two_alpha_R,[1 1 nbeams]);
ECHO = double(inadcp.ECHO(:,EnsIndex,:));                                  % echo intensity (counts)

SV = C + 10*log10((Tx+273.16).*R.^2) - LDBM - PDBW + two_alpha_R + Kc.*(ECHO-Er);
%SV = C + 10*log10((Tx+273.16).*R.^2) - LDBM - PDBW + two_alpha_R + 10*log10(10.^(Kc*ECHO/10)-10.^(Kc*Er/10))-20*log10(0.96/131);
