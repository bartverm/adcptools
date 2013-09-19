function isADCP=isADCPstruct(inp,varargin)
% ISADCPSTRUCT checks whether a variable is an ADCP structure
%       TF = isADCPstruct(VAR) Checks wheter variable VAR is an ADCP
%       structure.
%       VAR can be a structure, a structure array or a cell of structures.
%       in the latter two cases the function will return a boolean array
%       with the same size as the input
%
%       TF = isADCPstruct(VAR,FLAGS) speciefies in FLAGS on which data the
%       function should check. FLAGS is a character array with any of the
%       following characters (case-insensitive):
%       e   ensemble info
%       v   velocity data
%       c   correlation data
%       h   echo intensity
%       p   percentage good
%       b   bottom track data
%       x   external data stored in the adcp file
%       FLAGS is equal to 'ev' by default
%
%`      See Also: readadcp2
%
%       Author: Bart Vermeulen
%       Last edit: 05-08-2010

%    Copyright 2013 Bart Vermeulen
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

%% Input checking
if verLessThan('matlab','7.13')
    error(nargchk(1,2,nargin)) %#ok<NCHKN>
    error(nargoutchk(0,1,nargout)) %#ok<NCHKE>
else
    narginchk(1,2)
    nargoutchk(0,1)
end


if isstruct(inp) 
    isADCP=false(size(inp));
    inp={inp};
elseif iscell(inp)
    isADCP=false(size(inp));
else
    isADCP=false;
    return
end

if nargin==2 && ~ischar(varargin{1}), throw(MException('isADCPstruct:invalidInputFlagsClass','Flags should be given as a character array'))
else varargin{1}='ev';
end

%% Definition of fields to be present in adcp structure
fields={'firmver';'firmrev';'sysconf';'sysconfstr';'SymData';'LagLength';...
        'usedbeams';'nbins';'pingperens';'binsize';'blnk';'minthrsh';...
        'ncodrep';'minpercgood';'maxerrvel';'Tbetweenpng';'corinfo';...
        'corstr';'headalign';'headbias';'sensource';'senavail';...
        'distmidbin1';'lngthtranspulse';'watrefbins';'mintarget';...
        'lowlattrig';'distpulse';'cpuserial';'bandwidth';'syspower';...
        'basefreqid';'serial'};
if ~isempty(strfind(varargin{1},'e'))
fields=[fields; {'ensnum';'BITcheck';'speedsound';'depthtransd';'heading';'pitch';...
          'roll';'salinity';'temperature';'prepingT';'headstd';'pitchstd';...
          'rollstd';'ADC';'errorstat1';'errorstat2';'errorstat3';...
          'errorstat4';'pressure';'pressurevar';'timeV'}];
end
if ~isempty(strfind(varargin{1},'b'))
fields=[fields; {'btpingperens';'reacqdelay';'mincormag';'minevampl';...
          'btminpergood';'btmode';'btmaxerrv';'btrange';'btvel';'btcor';...
          'btevampl';'btpercgood';'minlyrsize';'nearbnd';'farbnd';...
          'reflyrvel';'reflyrcor';'reflyrint';'reflyrpergood';'maxdepth';...
          'rssiamp';'gain';'btrange'}];
end
if ~isempty(strfind(varargin{1},'v')),fields=[fields; {'VEL'}]; end
if ~isempty(strfind(varargin{1},'h')),fields=[fields; {'ECHO'}]; end
if ~isempty(strfind(varargin{1},'c')),fields=[fields; {'CORR'}]; end
if ~isempty(strfind(varargin{1},'p')),fields=[fields; {'PERC'}]; end
if ~isempty(strfind(varargin{1},'x')),fields=[fields; {'nfiles';'hfiles';'dfiles'}]; end

%% Check existance of the selected fields
for cc=1:numel(inp)
    if ~isstruct(inp{cc}), continue, end
    if isempty(strfind(varargin{1},'x')), isADCP(cc)=all(isfield(inp{cc},fields));
    else isADCP(cc)=all(isfield(inp{cc},fields(1:(end-3)))) && any(isfield(inp{cc},fields((end-2):end))); 
    end
end