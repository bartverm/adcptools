function misal=getExtMisalign(inadcp)
% misal=GetExtMisalign(INADCP)

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

if verLessThan('matlab','7.13')
    error(nargchk(1,1,nargin)) %#ok<NCHKN>
    error(nargoutchk(0,1,nargout)) %#ok<NCHKE>
else
    narginchk(1,1)
    nargoutchk(0,1)
end

misal=nan(size(inadcp));

if isstruct(inadcp), inadcp={inadcp(:)};end

for cc=1:numel(inadcp)
    if ~isADCPstruct(inadcp(cc)), continue, end
    extH=getADCPHeading(inadcp{cc});
    if isempty(extH), misal(cc)=nan; continue, end
    extV=[cosd(extH); sind(extH); extH*0]';
    intH=double(inadcp{cc}.heading(:))/100;
    if nanmax(intH)-nanmin(intH)<1e-6, misal(cc)=nan; continue, end
    intV=[cosd(intH),sind(intH),intH*0];
    dangle=atan2(sum(cross(extV,intV,2),2),dot(extV,intV,2));
    misal(cc)=atan2(nanmean(sin(dangle)),nanmean(cos(dangle)))/pi*180;
end
