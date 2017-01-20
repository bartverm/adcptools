function dat_out=readDeployment(searchstr)
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


allFiles=dir(searchstr); % read all filenames
assert(~isempty(allFiles),'readDeployment:NoFileFound','Could not find any file'); % check there is at least one file with the given deployment name
allFiles([allFiles.isdir])=[];


%% Read Transect Files
is_pd0_new=~cellfun(@isempty, regexp({allFiles.name},'.*\.PD0$'));
is_pd0_old=~cellfun(@isempty, regexp({allFiles.name},'.*[0-9]{3,3}r\.[0-9]{3,3}$'));
if any(is_pd0_new) && any(is_pd0_old)
    error('Found both old and new style filenames');
elseif ~any(is_pd0_new) && ~any(is_pd0_old)
    error('Could not find raw data files');
elseif any(is_pd0_new)
    is_pd0=is_pd0_new;
    is_new=true;
else
    is_pd0=is_pd0_old;
    is_new=false;
end
raw_files=cellfun(@fullfile, {allFiles(is_pd0).folder}, {allFiles(is_pd0).name},'UniformOutput', false);

dat_out=readADCP(raw_files); % Read files

%% Read Navigation Files
if is_new
    is_nav=~cellfun(@isempty, regexp({allFiles.name},'.*GPS\.TXT$'));
else
    is_nav=~cellfun(@isempty, regexp({allFiles.name},'.*[0-9]{3,3}n\.[0-9]{3,3}$'));
end
if any(is_nav)
    nav_files=cellfun(@fullfile, {allFiles(is_nav).folder}, {allFiles(is_nav).name},'UniformOutput', false);
end
nav_nmea=readNMEA(nav_files);
dat_out.nav_files=match_nmea_adcp(dat_out,nav_nmea);



% %% Read Depth Sounder Files
% dfiles=match_and_cat({'.*[0-9]{3,3}d\.[0-9]{3,3}$'; '.*DS\.TXT$'}); % search for depth sounder files
% if ~isempty(dfiles) % found something
%     disp('Reading depth sounder files...')
%     outADCP.dFiles=readNMEAADCP(outADCP,dfiles); % read it
% end
% 
% %% Read External Heading Files
% hfiles=match_and_cat({'.*[0-9]{3,3}h\.[0-9]{3,3}$'; '.*EH\.TXT$'}); % search for depth sounder files
% if ~isempty(hfiles)
%     disp('Reading external heading files...')
%     outADCP.hFiles=readNMEAADCP(outADCP,hfiles);
% end
% 
% %% Read Transect files
% tfiles=match_and_cat('.*[0-9]{3,3}t\.[0-9]{3,3}$'); % search for transect files
% if ~isempty(tfiles) 
%     disp('Reading transect files...')
%     outADCP.tFiles=readTfiles(outADCP,tfiles);
% end
% 
% %% Read VISEA Extern file
% vfiles=match_and_cat('.*[0-9]{3,3}extern\.dat$'); % search for transect files
% if ~isempty(vfiles)
%     disp('Reading VISEA extern files...')
%     outADCP.VISEA_Extern=readViseaExtern(outADCP,vfiles);
% end
% 
% %% Function to match regular expression in file names and concatenate the path to the result
% function files=match_and_cat(rex)
%     if ~iscellstr(rex) && ischar(rex) % Make a char input into 1x1 cell of char
%         rex={rex};
%     end
%     assert(iscellstr(rex)); % Make sure input is cell of char
%     files={}; % Initialize empty list of files
%     for crex=1:numel(rex) % Match each regular expression given
%         tmpfiles=regexp(allFiles,rex{crex},'match'); % search for cells matching regular expression
%         files=[files;tmpfiles{~cellfun(@isempty,tmpfiles)}]; %#ok<AGROW> % Get result as cell of strings and add to output(growth should be insignificant)
%     end
%     files=cellfun(@(x,y) [x y],repmat({path},size(files)),files,'UniformOutput',false); % Concatenate path and filenames
% end
% 
% end

