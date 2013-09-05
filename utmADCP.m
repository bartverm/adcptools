function [UTMx,UTMy,zone]=utmADCP(inadcp)
% Calculates UTM coordinates and zone for lat long coordinates in adcp
% files
%
% Requires mapping toolbox
% Requires external gps data in adcp file
% Assumes latitude and longitude are in the same zone

%    Copyright 2009,2010 Bart Vermeulen, Maximiliano Sassi
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

% Last edit 15-09-2009

% search for coordinates in adcp structure

gf1='';
gf2='';
UTMx=[];
UTMy=[];
zone='';

mfields={'nFiles','hFiles','dFiles'};
sfields={'GGA','RMC','GLL'};
cntmf=0;
cntsf=0;
while isempty(gf2)
    cntsf=cntsf+1;
    while isempty (gf1)
        cntmf=cntmf+1;
        if isfield(inadcp,mfields{cntmf}) && isfield(inadcp.(mfields{cntmf}),sfields{cntsf})
            gf1=mfields{cntmf};
            gf2=sfields{cntsf};
        end
        if cntmf==length(mfields)
            cntmf=0;
            break
        end
    end
    if cntsf==length(sfields)
        break;
    end
end
if ~isempty(gf1)
    lat=inadcp.(gf1).(gf2).lat;
    long=inadcp.(gf1).(gf2).long;
    if strcmp(gf2,'GGA')
        badGPS=inadcp.(gf1).(gf2).qualind~=1 & inadcp.(gf1).(gf2).qualind~=2 ;
        lat(badGPS)=NaN;
        long(badGPS)=NaN;
    else
        warning('utmADCP:QualInd','No Quality indicator found')
    end
else
    cntsf=0;
    while isempty (gf1)
        cntsf=cntsf+1;
        if isfield(inadcp,sfields{cntsf})
            gf1=sfields{cntsf};
        end
        if cntsf==length(sfields)
            break
        end
    end
    if ~isempty(gf1)
        lat=inadcp.(gf1).lat;
        long=inadcp.(gf1).long;
        if strcmp(gf1,'GGA')
            badGPS=inadcp.(gf1).qualind~=1 & inadcp.(gf1).qualind~=2 ;
            lat(badGPS)=NaN;
            long(badGPS)=NaN;
        else
            warning('utmADCP:QualInd','No Quality indicator found')
        end
    else
        if isfield(inadcp,'tFiles') && all(isfield(inadcp.tFiles,{'lat','long'}))
            lat=inadcp.tFiles.lat;
            long=inadcp.tFiles.long;
        else
            error('utmADCP:NoData','Could not find coordinates information in ADCP file')
        end
    end
end


% find utm zone for given data
goodff=(~isnan(lat)) & (~isnan(long));
if ~any(goodff), return; end
zone=utmzone(lat(goodff),long(goodff));

% create map structure
ellipsoid = utmgeoid(zone);
adcpmap = defaultm('utm'); 
adcpmap.zone = zone; 
adcpmap.geoid = ellipsoid; 
adcpmap.flatlimit = []; 
adcpmap.maplatlimit = []; 
adcpmap = defaultm(adcpmap);
[UTMx,UTMy]=mfwdtran(adcpmap,lat,long);
UTMx=UTMx(:);
UTMy=UTMy(:);
