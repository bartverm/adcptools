function heading=getADCPHeading(inadcp)
% heading=getADCPHeading(INADCP)

%    Copyright 2010 Bart Vermeulen
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

gf1='';
gf2='';

mfields={'nFiles','hFiles','dFiles'};
sfields={'HPR','HDT'};
cntmf=0;
cntsf=0;
if isfield(inadcp,'hFiles')
    warning('getADCPHeading:HeadInRFile','There is data from external ascii heading files, \n so external heading is recorded in the heading field of the data structure')
end
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
    heading=inadcp.(gf1).(gf2).heading;
elseif isfield(inadcp,'VISEA_Extern') && isfield(inadcp.VISEA_Extern,'Heading')
    heading=inadcp.VISEA_Extern.Heading;
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
        heading=inadcp.(gf1).heading;
    else
        heading=double(inadcp.heading)/100;
    end
end
heading=heading(:)';

