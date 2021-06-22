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
%     if strcmp(gf2,'GGA')
%         badGPS=inadcp.(gf1).(gf2).qualind~=1 & inadcp.(gf1).(gf2).qualind~=2 ;
%         lat(badGPS)=NaN;
%         long(badGPS)=NaN;
%     else
%         warning('utmADCP:QualInd','No Quality indicator found')
%     end
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
            warning('utmADCP:NoGPS','Could not find coordinates information in ADCP file, using bottom tracking')
            if ~isfield(inadcp,'btvel')
                error('utmADCP:NoGPSNoBt','Could not find bottom tracking data')
            end
            [inadcp.VEL,inadcp.btvel] = filterADCP(inadcp,'','filterBT',true); % Filter velocity
            [~, btvel]=corADCP(inadcp,'e'); % transform to earth velocity
            btvel=(btvel(1:end-1,1:2)+btvel(2:end,1:2))/2; % average two consecutive bottom tracking velocity (matches velocity with delta time)
            time=datenum(inadcp.timeV); % Get the time
            
            % STREAMPRO hack for time vector
            if ~all(diff(time)>0) % If time is not monotonicallyincreasing
                tv=inadcp.timeV1C; % Use non 2K compliant time
                tv(:,1)=tv(:,1)+2000; % add century (assuming it is 21st century)
                time=datenum(tv); % store result in time array
            end
            % END STREAMPRO hack for time vector
            
            dt=diff(time)*24*3600; % compute delta time
            fgood=all(isfinite(btvel),2); % find bad bottom tracking
            btvel(~fgood,:)=interp1(time(fgood),btvel(fgood,:),time(~fgood),'pchip'); % interpolate missing bottom tracking
            pos=cumsum([0 0; bsxfun(@times,btvel,dt)],1); % compute position based on bottom tracking
            
            % Reset position to (0,0) at the start of a new file (avoids
            % strange jumps in position due to large delta_t)
            f_newfile=find(diff(inadcp.FileNumber)>0); % find where files end
            reset_pos=zeros(size(pos)); % vector to hold start position of each track
            find_reset_pos=inadcp.FileNumber'-1; % vector to hold indices with file number for each ensemble
            reset_pos(find_reset_pos>0,:)=pos(f_newfile(find_reset_pos(find_reset_pos>0))+1,:); % Store for each file the starting coordinates
            pos=pos-reset_pos; % Subtract starting coordinates for each file
            
            % store output and return
            UTMx=pos(:,1);
            UTMy=pos(:,2);
            return
        end
    end
end

[UTMx UTMy zone] = geo2utm(lat,long);

UTMx=UTMx(:);
UTMy=UTMy(:);

