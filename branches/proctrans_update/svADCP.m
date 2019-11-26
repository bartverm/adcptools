function [SV,Kc,C] = svADCP(inadcp,pt,varargin)
% svADCP compute signal backscatter strength (dB)
%   SV = svADCP(ADCP, PT) computes signal backscatter strength from 
%       ADCP profile data. ADCP is an ADCP structure as read by 
%       readADCP. PT is an acoustics.PistonTransducer object which
%       specifies the transducer and water conditions. This can be
%       generated manually or from the ADCP structure with the
%       piston_transducer_from_adcpstruct function. This function
%       computes backscatter strength ignoring sediment attenuation and
%       uses the sonar equation as proposed by Gostiaux and van Haren (
%       2010) which is an improved version of the sonar equation by 
%       Deines (1999)
%
%   SV = svADCP(ADCP, PT, EnsRange) compute backscatter strenght for
%       selected ensemble range, EnsRange is a two column vector with
%       starting and ending indexes
%
%   SV = svADCP(...,'ParameterName', ParameterValue) allows to specify the
%       following optional settings:
%       
%       'Kc'
%       Scalar value indicating the RSSI scale factor dB per count
%       (default: 0.43). This value can be calibrated with a hydrophone or
%       see Sassi et al. 2012 for an alternative method
%
%       'C' 
%       Scalar value indicating the instrument constant. If not specified C
%       is estimated based on instrument characteristics (Deines, 1999)
%
%       'Er'
%       RSSI in still water (noise level). Default is 40. See PT3 command
%
%       'T_offset'
%       Temperature sensor offset, see PS0 command. Default is -0.35
%       
%  See also: readadcp, piston_transducer_from_adcpstruct
%

%    Copyright 2007-2019 Bart Vermeulen, Maximiliano Sassi
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


%% Parsing input
P = inputParser;                                                           % Create an input parser
P.FunctionName = 'svADCP';                                                 % Check for correct spelling in function
P.addRequired('inadcp',@isstruct);                                         % Check adcp structure
P.addRequired('pt',@(x) isscalar(x) && isa(x,'acoustics.PistonTransducer'))% Piston Transducer properties
P.addOptional('EnsRange',[0 0],@(x) isnumeric(x) && numel(x)==2);          % Ensemble range to process
P.addParameter('Kc',0.43,@(x) isscalar(x) && isnumeric(x) && x>0);         % RSSI scale factor dB per count (calibration !) Note that Aquavision uses 0.43
P.addParameter('Er',40,@(x) isscalar(x) && isnumeric(x) && x>0);           % RSSI in a bucket, see PT3 command
P.addParameter('T_Offset',-0.35,@(x) isscalar(x) && isnumeric(x));         % Temp Sens Offset (PS0 results)
P.addParameter('C',0,@(x) isscalar(x) && isnumeric(x));

P.parse(inadcp,pt,varargin{:});                                               % Parse input

EnsRange = P.Results.EnsRange;
% IsHadcp = P.Results.IsHadcp;
Kc= P.Results.Kc;                                                    
Er = P.Results.Er;                                                                   
t_offset = P.Results.T_Offset;                         
C=P.Results.C; 

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
% constants from rdi manuals
current_fact=11451/1e6;
switch inadcp.sysconf(1,1:3) 
    case '000'
        volt_fact=2092719/1e6;
        current_fact=43838/1e6;
        Ctmp=-159.1;
    case '100'
        volt_fact=592157/1e6;
        Ctmp=-153.0;
    case '010'
        volt_fact=592157/1e6;
        Ctmp=-143;
    case '110'
        volt_fact=380667/1e6;
        Ctmp=-139.3;
    case '001'
        volt_fact=253765/1e6;
        Ctmp=-129.1;
    case '101'
        volt_fact=253765/1e6;
        Ctmp=NaN;
end
if any(strcmp(P.UsingDefaults,'C'))
    C=Ctmp;
end

%% Step 3: obtain selected inadcp parameters
B = double(inadcp.blnk(inadcp.FileNumber(EnsIndex)))/100;                  % blank after transmit (m)
L = double(inadcp.lngthtranspulse(inadcp.FileNumber(EnsIndex)))/100;       % transmit pulse length (m)
LDBM = 10*log10(L);                                                         
D = double(inadcp.binsize(inadcp.FileNumber(EnsIndex)))/100;               % depth cell length (m)
curr = reshape(current_fact*double(inadcp.ADC(EnsIndex,1)),1,[]);          % current (A)
volt = reshape(volt_fact*double(inadcp.ADC(EnsIndex,2)),1,[]);             % voltage (V)
DC_COEF = 9.82697464e1;                                                    % Temperature coefficients
FIRST_COEF = -5.86074151382e-3;
SECOND_COEF = 1.60433886495e-7;
THIRD_COEF = -2.32924716883e-12;
t_cnts = reshape(double(inadcp.ADC(EnsIndex,6))*256,1,[]);                 % Temperature Counts (ADC value)
Tx = t_offset + ((THIRD_COEF.*t_cnts + SECOND_COEF).*t_cnts +...
    FIRST_COEF).*t_cnts + DC_COEF;                                         % real-time temperature of the transducer (C)

%% Step 4: prescribe relevant external variables
%TODO: deal with stratified conditions (variable speed of sound)

%% Step 5: determine transmit power
power = volt.*curr;                                                        % Transmit power in Watts
PDBW = 10*log10(power);

%% Step 6: evaluate variables in sonar equation for each depth cell
R = (repmat(B + (D+L)/2,nbins,1) + repmat((1:nbins)'-1,1,nens).*...
    repmat(D,nbins,1) + repmat(D,nbins,1)/4)./repmat(cosd(bangle),nbins,1);%*creal/c  if correction for speed of sound is needed
                                                                           % creal is the average sound speed from transducer to the range cell.
alpha=pt.attenuation;
alpha_n = 2*alpha*D./cosd(bangle);                                         % absorption for each range cell
two_alpha_R = 2*alpha*repmat(B,nbins,1)./repmat(cosd(bangle),nbins,1) +...
    repmat((1:nbins)',1,nens).*repmat(alpha_n,nbins,1);                    % compute 2alphaR

%% Step 7: calculate profiles of backscatter coefficient
ECHO = double(inadcp.ECHO(:,EnsIndex,:));                                  % echo intensity (counts)
SV = C + 10*log10((Tx+273.16).*R.^2.*pt.near_field_correction(R).^2) - LDBM - PDBW + two_alpha_R + 10*log10(10.^(Kc.*ECHO/10)-10.^(Kc.*Er/10));
