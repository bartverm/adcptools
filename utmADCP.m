function [UTMx,UTMy,zone]=utmADCP(inadcp)
% Calculates UTM coordinates and zone for lat long coordinates in adcp
% files
%   
%   [UTMx,UTMy,zone]=utmADCP(inadcp) converts lat long coordinates (assumed
%       to be measured on WGS84 into UTMx and UTMy coordinates, in the
%       projected UTM coordinate system. Also the zone is returned based on
%       the average location in the data. 
%

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
        if isfield(inadcp,'VISEA_Extern') && all(isfield(inadcp.VISEA_Extern,{'Latitudeseconds','Longitudeseconds'}))
            lat=inadcp.VISEA_Extern.Latitudeseconds'/3600;
            long=inadcp.VISEA_Extern.Longitudeseconds'/3600;
        elseif isfield(inadcp,'tFiles') && all(isfield(inadcp.tFiles,{'lat','long'}))
            lat=inadcp.tFiles.lat';
            long=inadcp.tFiles.long';
        else
            error('utmADCP:NoData','Could not find coordinates information in ADCP file')
        end
    end
end

%% find utm zone for given data
% select good data
goodff=(~isnan(lat)) & (~isnan(long));
if ~any(goodff), return; end

% find zone
mlat=nanmean(lat);
mlong=nanmean(long);
lts = [-80:8:72 84]';
lns = (-180:6:180)';
latzones = char([67:72 74:78 80:88]');

indx = find(lts <= mlat);
ltindx = indx(max(indx));

indx = find(lns <= mlong);
lnindx = indx(max(indx));

if ltindx < 1 || ltindx > 21
    ltindx = [];
elseif ltindx == 21
    ltindx = 20;
end

if lnindx < 1 || lnindx > 61
    lnindx = [];
elseif lnindx == 61
    lnindx = 60;
end

zone = [num2str(lnindx) latzones(ltindx)];



% transform with suitable package
% if exist('utm_fwd','file')==2 % GeographicLib is available
%     disp('Using GeographicLib for UTM forward transformation')
%     if mean(lat(goodff))>0, northp=true; else northp=false; end
%     [UTMx,UTMy]=utm_fwd(lnindx,northp,lat,long);
% elseif license('checkout','map_toolbox')
%     disp('Using Mapping Toolbox for UTM forward transformation')
% %     zone=utmzone(lat(goodff),long(goodff));
%     ellipsoid = utmgeoid(zone);
%     adcpmap = defaultm('utm'); 
%     adcpmap.zone = zone; 
%     adcpmap.geoid = ellipsoid; 
%     adcpmap.flatlimit = []; 
%     adcpmap.maplatlimit = []; 
%     adcpmap = defaultm(adcpmap);
%     [UTMx,UTMy]=mfwdtran(adcpmap,lat,long);
% else
%     disp('Using custom UTM forward transformation')
    lns = (-180:6:180)';
    indx = find(lns <= mean(long(goodff)));
    zone = indx(max(indx));
    if zone < 1 || zone > 61
        zone = [];
    elseif zone == 61
        zone = 60;
    end
    long0=(zone*6-183)/180*pi;
    long=long/180*pi;
    lat=lat/180*pi;
    a=6378.137;
    finv=298.257223563; f=1/finv;
    if mean(lat(goodff))>0, N0=0; else N0=10000; end
    k0=0.9996;
    E0=500;
    n=f/(2-f);
    t=sinh(atanh(sin(lat))-2*sqrt(n)/(1+n)*atanh(2*sqrt(n)*sin(lat)/(1+n)));
    zeta=atan(t./cos(long-long0));
    eta=atanh(sin(long-long0)./sqrt(1+t.^2));
    A=a/(1+n)*(1+n.^2/4+n^4/64);
    alph(1)=1/2*n-2/3*n^2+5/16*n^3;
    alph(2)=13/48*n^2-3/5*n^3;
    alph(3)=61/240*n^3;
    jvec=1:3;
    UTMx=E0+k0*A*(eta+sum(bsxfun(@times,alph,cos(2*bsxfun(@times,jvec,zeta)).*sinh(2*bsxfun(@times,jvec,eta))),2));
    UTMy=N0+k0*A*(zeta+sum(bsxfun(@times,alph,sin(2*bsxfun(@times,jvec,zeta)).*cosh(2*bsxfun(@times,jvec,eta))),2));  
    UTMx=UTMx*1000;
    UTMy=UTMy*1000;
% end

UTMx=UTMx(:);
UTMy=UTMy(:);
