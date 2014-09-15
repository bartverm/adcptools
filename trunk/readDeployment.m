function [outADCP]=readDeployment(DepName,varargin)
% READDEPLOYMENT Read a Winriver deployment
%   ADCP = readDeployment(DEPNAME) read deployment files starting with
%       DEPNAME in the current folder. It currently supports navigation,
%       depth-sounder, external heading and transect files.
%   
%   ADCP = readDeployment(DEPNAME,PATH) reads files in PATH
%
% See also readADCP

%    Copyright 2009,2014 Bart Vermeulen
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


%% parse input 

%path
if nargin>1 % if path is given
    path=varargin{1}; % get it
    assert(ischar(path) && isrow(path),'readDeployment:BadFormatPath','Path must be a row vector of characters') % make sure it is a row vector of characters
    assert(exist(path,'dir')==7,'readDeployment:InexistentPath',['Cannot find folder: ',path]) % check the path exists

    % append slash to path if it's missing
    if ~any(strcmp(path(end),{'/','\'})) % is a trailing slash or backslash missing?
        if ispc(), path=[path,'\']; else path=[path,'/'];end % if we are on windows add a backslash, otherwise add a slash
    end % if
else % if path is not given 
    path=''; % make path empty char
end % if

% Depname
allFiles=dir([path,DepName,'*']); % read all filenames
assert(~isempty(allFiles),'readDeployment:NoFileFound',['Could not find any file for deployment: ', DepName]); % check there is at least one file with the given deployment name
allFiles=({allFiles(~[allFiles(:).isdir]).name})';


%% Read Transect Files
rfiles=match_and_cat({'\w*[0-9]{3,3}r\.[0-9]{3,3}'; '\w*.PD0'});% search for raw data files
assert(~isempty(rfiles),'ReadDeployment:NoRFiles','Could not find raw data files') % Make sure we found at least one adcp file
disp('Reading ADCP raw data files...') % Tell something nice to the user
outADCP=readADCP(rfiles); % Read files

%% Read Navigation Files
nfiles=match_and_cat({'\w*[0-9]{3,3}n\.[0-9]{3,3}'; '\w*GPS.TXT'}); % search for navigation files
if ~isempty(nfiles) % Found something?
    disp('Reading navigation files...') % Tell something nice to the user
    outADCP.nFiles=readNMEAADCP(outADCP,nfiles); % try to read it
end

%% Read Depth Sounder Files
dfiles=match_and_cat({'\w*[0-9]{3,3}d\.[0-9]{3,3}'; '\w*DS.TXT'}); % search for depth sounder files
if ~isempty(dfiles) % found something
    disp('Reading depth sounder files...')
    outADCP.dFiles=readNMEAADCP(outADCP,dfiles); % read it
end

%% Read External Heading Files
hfiles=match_and_cat({'\w*[0-9]{3,3}h\.[0-9]{3,3}'; '\w*EH.TXT'}); % search for depth sounder files
if ~isempty(hfiles)
    disp('Reading external heading files...')
    outADCP.hFiles=readNMEAADCP(outADCP,hfiles);
end

%% Read Transect files
tfiles=match_and_cat('\w*[0-9]{3,3}t\.[0-9]{3,3}'); % search for transect files
if ~isempty(tfiles) 
    disp('Reading transect files...')
    outADCP.tFiles=readTfiles(outADCP,tfiles);
end

%% Read VISEA Extern file
vfiles=match_and_cat('\w*[0-9]{3,3}extern\.dat'); % search for transect files
if ~isempty(vfiles)
    disp('Reading VISEA extern files...')
    outADCP.VISEA_Extern=readViseaExtern(outADCP,vfiles);
end

%% Function to match regular expression in file names and concatenate the path to the result
function files=match_and_cat(rex)
    if ~iscellstr(rex) && ischar(rex) % Make a char input into 1x1 cell of char
        rex={rex};
    end
    assert(iscellstr(rex)); % Make sure input is cell of char
    files={}; % Initialize empty list of files
    for crex=1:numel(rex) % Match each regular expression given
        tmpfiles=regexp(allFiles,rex{crex},'match'); % search for cells matching regular expression
        files=[files;tmpfiles{~cellfun(@isempty,tmpfiles)}]; %#ok<AGROW> % Get result as cell of strings and add to output(growth should be insignificant)
    end
    files=cellfun(@(x,y) [x y],repmat({path},size(files)),files,'UniformOutput',false); % Concatenate path and filenames
end

end
