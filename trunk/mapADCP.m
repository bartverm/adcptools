function [outx,outy,outz]=mapADCP(inadcp,varargin)
% MAPADCP determines the offset of the ensonified regions by the (H)ADCP 
%         with respect to the (H)ADCP centre. Definitions of X,Y,Z are
%         based on the RDI manuals.
%         [DX,DY,DZ]=mapADCP(INADCP) determines the offset of the ensonified
%         areas in x, y and z direction respectively East, North and Up
%         
%         [DX,DY,DZ]=mapADCP(INADCP,'PropertyName','PropertyValue') allows
%         to specify the following additional settings:
%
%         IsHADCP 
%         true | {false}
%         Logical value that sets if the structure contains data of a HADCP
%
%         bAngle
%         Numerical value that sets the beam-angle in degrees. If it's not
%         set the script will try to read it from the file. If that fails
%         too it will generate an error
%
%         avgHead
%         true | {false}
%         If set to true the function will output the offsets with respect
%         to the average heading of the instrument. This setting only works
%         if IsHADCP is set to true.
%
%         IsUpward
%         {true} | false
%         Determines wheter the instrument looks upward or downward. For
%         the HADCP upward corresponds with an upward oriented pressure
%         sensor
%
%         See also: readadcp2 depthADCP utmADCP
%

%    Copyright 2009-2010 Bart Vermeulen, Maximiliano Sassi
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

         
%         Last edit:  12-10-2010, allow pitch and roll to be as high as 40
%                     degrees

%
%         Last edit:  17-08-2010, calculation of pitch from raw tilts moved
%                     such that it's calculated before 180 degrees is added
%                     to the roll for upward looking case

%         Last edit: 17-05-2010 for HADCP beams are defined in downward
%         looking case (pressure sensor down) and afterwards rotated for
%         upward looking. 

%         Last edit: 10-07-2009 some minor corrections in cartesian
%         coordinates transformation with HADCP

%         Last edit: 08-07-2009 compute coordinates using maximum number of
%         bins, clean some memory while doing operations

%         Last edit: 06-07-2009


%% Handling input
P=inputParser;                                                             %Creating an input parser
P.addRequired('inadcp',@isstruct)
P.addParamValue('IsHADCP',false,@(x) islogical(x) && isscalar(x))
P.addParamValue('bAngle',NaN,@(x) isnumeric(x) && isscalar(x))
P.addParamValue('avgHead',false,@(x) islogical(x) && isscalar(x))
P.addParamValue('IsUpward',true,@(x) islogical(x) && isscalar(x))
P.addParamValue('UseExtHeading',false,@(x) islogical(x) && isscalar(x))
P.addParamValue('Beam3misalign',double(inadcp.headalign(inadcp.FileNumber))/100,@isnumeric)         % Misalignment in degrees (wind directions, only necessary with ext. heading)

P.parse(inadcp,varargin{:})

inadcp=P.Results.inadcp;
IsHADCP=P.Results.IsHADCP;
bangle=P.Results.bAngle/180*pi;
avgHead=P.Results.avgHead;
IsUpward=P.Results.IsUpward;
extHead=P.Results.UseExtHeading;
Beam3mis=P.Results.Beam3misalign;

clear P                                                                    % clear input parse

%% Determine beam angle
if isnan(bangle) 
    if IsHADCP
        bangle=max(double(inadcp.HADCPbeamangle)/180*pi);
        if bangle==0
            warning('mapADCP:UnknownBangle','Unknown beam-angle, Please provide one!')
        end
    else
        switch inadcp.sysconf(9:10)
            case '00'
            bangle=15/180*pi;
            case '10'
            bangle=20/180*pi;
            case '11'
            bangle=30/180*pi;
            case '01'
                error('mapADCP:UnknownBangle','Unknown beam-angle, Please provide one');
        end
    end
end

%% Change angles to radians and remove invalid angles
if extHead && ~IsHADCP
    heading=getADCPHeading(inadcp)';
    if isempty(heading)
        heading=double(inadcp.heading)/100;                                    %Read heading from internal sensor
    end
    heading=heading+Beam3mis;
else
    heading=double(inadcp.heading)/100;                                        %Change heading to doubles
end
heading(heading>359.99)=NaN;                                               %Filter out invalid headings
heading=heading/180*pi;                                                    %Change to radians
% the heading average option will rotate the coordinate system in such a
% way that the correspondence (x,y,z) -> (E,N,U) will be (x,y,z) -> (N,E,U)
if avgHead && IsHADCP
    mheading=atan2(nanmean(sin(heading)),nanmean(cos(heading)));
    heading=heading-mheading; % check !
end
clear mheading
pitch=double(inadcp.pitch)/100/180*pi;                                     %Change pitch to radians
pitch(abs(pitch)>2*pi/9)=NaN;                                                %Filter out invalid Pitches
roll=double(inadcp.roll)/100/180*pi;                                       %Change pitch to radians
roll(abs(roll)>2*pi/9)=NaN;                                                  %Filter out invalid Pitches
pitch=atan(tan(pitch).*cos(roll));                                         %Calculate real pitch from raw pitch (see adcp cor transf manual)
if IsUpward                                                                
    roll=roll+pi;
end

%% Determine Range of the instrument (z component)
distmidbin1=double(inadcp.distmidbin1)/100;                                %Get the distance of the middle of the first bin
binsize=double(inadcp.binsize)/100;                                        %Get the binsize
nbins=double(inadcp.nbins);                                                %Get the number of bins
nens=length(inadcp.ensnum);                                                %Find number of ensembles
FileNumber=repmat(inadcp.FileNumber,max(nbins),1);                         %Replicate the Number of dataset to the bins
BinDist=distmidbin1(FileNumber)+repmat((0:max(nbins)-1)',1,nens).*...
    binsize(FileNumber);                                                   %Calculate distance to each measured cell
nbins = max(nbins);                                                        % take the maximum number of bins
clear inadcp                                                               % clear memory

%% Calculate cartesian coordinates for untilted case
if IsHADCP
    yy=cat(3,repmat(BinDist,[1,1,3]),zeros(nbins,nens,1));                 % Y-axis is positive away from the HADCP
    tbangle=tan(bangle);                                                   % Beams are defined in the downward looking case and afterwards rotated for upward case
    xx(:,:,1)=-yy(:,:,1)*tbangle;
    xx(:,:,2)=yy(:,:,2)*tbangle;
    xx(:,:,3:4)=zeros(nbins,nens,2);
    zz=yy*0;
else
    zz=repmat(-BinDist,[1,1,4]);
    tbangle=tan(bangle);
    xx(:,:,1)=zz(:,:,1)*tbangle;
    xx(:,:,2)=-zz(:,:,2)*tbangle;
    xx(:,:,3:4)=zeros(nbins,nens,2);
    yy(:,:,1:2)=zeros(nbins,nens,2);
    yy(:,:,3)=-zz(:,:,3)*tbangle;
    yy(:,:,4)=zz(:,:,4)*tbangle;
end

clear BinDist FileNumber                                                   % clear memory


%% Apply the tilts 
% Find rotation matrix (see adcp coordinate transformation manual)
ch=cos(heading);
sh=sin(heading);
cp=cos(pitch);
sp=sin(pitch);
cr=cos(roll);
sr=sin(roll);

M11=ch.*cr+sh.*sp.*sr;                                                     % compute rotation matrix
M12=sh.*cp;
M13=ch.*sr-sh.*sp.*cr;
M21=-sh.*cr+ch.*sp.*sr; 
M22=ch.*cp;
M23=-sh.*sr-ch.*sp.*cr;
M31=-cp.*sr;
M32=sp;
M33=cp.*cr;

clear s* c* heading pitch roll                                             % clear memory

xx=squeeze([xx(:,:,1);xx(:,:,2);xx(:,:,3);xx(:,:,4)]);
yy=squeeze([yy(:,:,1);yy(:,:,2);yy(:,:,3);yy(:,:,4)]);
zz=squeeze([zz(:,:,1);zz(:,:,2);zz(:,:,3);zz(:,:,4)]);

% replicate rotation matrix while doing the product saves some memory
xxt = xx.*repmat(M11,[4*nbins,1])+yy.*repmat(M12,[4*nbins,1])+zz.*repmat(M13,[4*nbins,1]);
yyt = xx.*repmat(M21,[4*nbins,1])+yy.*repmat(M22,[4*nbins,1])+zz.*repmat(M23,[4*nbins,1]);
zzt = xx.*repmat(M31,[4*nbins,1])+yy.*repmat(M32,[4*nbins,1])+zz.*repmat(M33,[4*nbins,1]);

clear xx yy zz

%% Output result
outx=cat(3,xxt(0*nbins+(1:nbins),:),xxt(1*nbins+(1:nbins),:),xxt(2*nbins+(1:nbins),:));
outy=cat(3,yyt(0*nbins+(1:nbins),:),yyt(1*nbins+(1:nbins),:),yyt(2*nbins+(1:nbins),:));
outz=cat(3,zzt(0*nbins+(1:nbins),:),zzt(1*nbins+(1:nbins),:),zzt(2*nbins+(1:nbins),:));
if ~IsHADCP
    outx=cat(3,outx,xxt(3*nbins+(1:nbins),:));
    outy=cat(3,outy,yyt(3*nbins+(1:nbins),:));
    outz=cat(3,outz,zzt(3*nbins+(1:nbins),:));
end
