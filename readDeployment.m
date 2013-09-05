function [outADCP]=readDeployment(DepName,path,varargin)
% READDEPLOYMENT Read a Winriver deployment
%                ADCP = readDeployment(DEPNAME,PATH) read deployment files
%                DEPNAME contained in PATH folder. Search and read
%                navigation, depth sounder and external heading files.
%                Works for WinRiver and WinRiver II files.
%
% See also readADCP

%    Copyright 2009 Bart Vermeulen
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

%% Read Transect Files
rfiles=dir([path,DepName,'*r.*']);                                         % Look for WinRiver transect files
rfiles2 = dir([path,DepName,'*.PD0']);                                     % Look for WinRiverII transect files
if isempty(rfiles)&&isempty(rfiles2)                                       % Stop running
    error('ReadDeployment:NoRFiles','Could not find raw data files')
end
disp('Reading ADCP raw data files...')
if isempty(rfiles)&&~isempty(rfiles2)                                      % Choose the ~empty structure
    rfiles=vertcat(rfiles2(:).name);
elseif isempty(rfiles2)&&~isempty(rfiles)
    rfiles=vertcat(rfiles(:).name);
end
rfiles= [repmat(path,[size(rfiles,1),1]),rfiles];
outADCP=readADCP(rfiles,varargin{:});
clear rfiles rfiles2

%% Read Navigation Files
nfiles=dir([path,DepName,'*n.*']);                                         % Look for WinRiver navigation files
nfiles2 = dir([path,DepName,'*GPS.TXT']);                                  % Look for WinRiverII navigation files
if isempty(nfiles)&&isempty(nfiles2)                                       % Continue running
    warning('ReadDeployment:NoNFiles','Could not find external navigation files')
else
    disp('Reading WinRiver navigation files...')
    if isempty(nfiles)&&~isempty(nfiles2)                                  % Choose the ~empty structure
        nfiles=vertcat(nfiles2(:).name);
    elseif isempty(nfiles2)&&~isempty(nfiles)
        nfiles=vertcat(nfiles(:).name);
    end    
    nfiles= [repmat(path,[size(nfiles,1),1]),nfiles];
    outADCP.nFiles=readNMEAADCP(outADCP,nfiles);
end
clear nfiles nfiles2

%% Read Depth Sounder Files
dfiles=dir([path,DepName,'*d.*']);                                         % Look for WinRiver depth sounder files
dfiles2 = dir([path,DepName,'*DS.TXT']);                                   % Look for WinRiverII depth sounder files
if isempty(dfiles)&&isempty(dfiles2)                                       % Continue running
    warning('ReadDeployment:NoDFiles','Could not find external depth files')
else
    disp('Reading WinRiver depth files...')
    if isempty(dfiles)&&~isempty(dfiles2)                                  % Choose the ~empty structure
        dfiles=vertcat(dfiles2(:).name);
    elseif isempty(dfiles2)&&~isempty(dfiles)
        dfiles=vertcat(dfiles(:).name);
    end    
    dfiles= [repmat(path,[size(dfiles,1),1]),dfiles];
    outADCP.dFiles=readNMEAADCP(outADCP,dfiles);
end
clear dfiles dfiles2

%% Read External Heading Files
hfiles=dir([path,DepName,'*h.*']);                                         % Look for WinRiver external heading files
hfiles2 = dir([path,DepName,'*EH.TXT']);                                   % Look for WinRiverII external heading files
if isempty(hfiles)&&isempty(hfiles2)                                       % Continue running
    warning('ReadDeployment:NoHFiles','Could not find external heading files')
else
    disp('Reading WinRiver heading files...')
    if isempty(hfiles)&&~isempty(hfiles2)                                  % Choose the ~empty structure
        hfiles=vertcat(hfiles2(:).name);
    elseif isempty(hfiles2)&&~isempty(hfiles)
        hfiles=vertcat(hfiles(:).name);
    end
    hfiles= [repmat(path,[size(hfiles,1),1]),hfiles];
    outADCP.hFiles=readNMEAADCP(outADCP,hfiles);
end

%% Read Transect files
tfiles=dir([path,DepName,'*t.*']);                                         % Look for WinRiver external heading files
% hfiles2 = dir([path,DepName,'*EH.TXT']);                                   % Look for WinRiverII external heading files
if isempty(tfiles) %&&isempty(hfiles2)                                       % Continue running
    warning('ReadDeployment:NoTFiles','Could not find transect files')
else
    disp('Reading WinRiver transect files...')
%     if isempty(tfiles)&&~isempty(hfiles2)                                  % Choose the ~empty structure
%         hfiles=vertcat(hfiles2(:).name);
%     elseif isempty(hfiles2)&&~isempty(hfiles)
     tfiles=vertcat(tfiles(:).name);
%     end
    tfiles= [repmat(path,[size(tfiles,1),1]),tfiles];
    outADCP.tFiles=readTfiles(outADCP,tfiles);
end
clear tfiles
