function out=readViseaExtern(adcp,filenames,varargin)
% READVISEAEXTERN Reads Visea extern files
%   READVISEAEXTERN(ADCP,FILENAME) Reads the visea file given by FILENAME
%       and returns a structure with the data stored in that file. It will
%       match the visea data and the ADCP data based on the ensemble
%       numbers in the file and the ADCP structure
%
%   READVISEAEXTERN(...,'ParamName',ParamValue) Allows to set the
%       following parameters:
%
%       'FieldWidth'
%       A scalar positive integer indicating the field width in the visea
%       files (default is 25)
%
%   See also: readADCP 

%    Copyright 2014 Bart Vermeulen
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
inp=inputParser;     % Create an object of the InputParser class
inp.addRequired('adcp',@isstruct);                                       % Add required input inadcp
inp.addRequired('filenames',@(x) (iscellstr(x) | ischar(x)));           % Add the required variable 'mmeafilename' and check for right format
inp.addParamValue('FieldWidth',25,@(x) isnumeric(x) && isscalar(x) && x>0 && floor(x)==x);    % Param value to set the field width in the visea file (it looks always like 25 chars wide)
inp.parse(adcp,filenames,varargin{:});                                            % Parse input
if ischar(filenames)
    filenames=cellstr(filenames);                                    % Change character array to cell
end
fwidth=inp.Results.FieldWidth;
clear inp

%% read data
% Initialize some variables
nFiles=numel(filenames); % Number of files to be read
out=struct(); % Initialize the output structure
out.timeV=nan(size(adcp.timeV)); % Add time vector field

for cf=1:nFiles
    fid=fopen(filenames{cf},'r'); % open file
    assert(fid~=-1,'readViseaExtern:FileNotOpen',['Could not open file:',filenames{cf}]); % check file has opened successfully
    head=fgetl(fid); % read header line
    assert(strcmp(head(1:42),'[YYYY/MM/DD, HH:MM:SS.SSS] ENSEMBLE NUMBER'),'readViseaExtern:BadHead',['Could not read Visea exter file: ',filenames{cf}]) % Check the header looks familiar
    head(1:43)=[]; % Remove fixed part of header
    head(end)=[]; % Remove newline character
    nFields=round(numel(head)/fwidth); % estimate number of fields
    head=mat2cell(reshape(head(1:min(end,nFields*fwidth)),[],nFields)',ones(nFields,1),fwidth); % Get field names from header
    head=cellfun(@(x) x(isstrprop(x,'alphanum')),head,'UniformOutput',false); % Retain only alphanumeric characters of field names (to make them suitable as variable name)
    dat=textscan(fid,['%*1s%4f/%2f/%2f,%2f:%2f:%6f%*1sENSEMBLE%5f:',repmat(['%',num2str(fwidth),'f'],1,nFields)]); % Read in the data
    fclose(fid); % close file
    
    % match ensnum
    [~,idx_visea,idx_adcp]=intersect(dat{7},adcp.ensnum); % find matching ensemble numbers
    tmptime=[dat{1:6}]; % Temporarily store time from Visea file
    out.timeV(idx_adcp,:)=tmptime(idx_visea,:); % Store visea time in output structure

    for cfld=1:nFields % loop over all fields
        if ~isfield(out,head{cfld}) % If field is not yet in output structure, add it 
           out.(head{cfld})=nan(1,size(out.timeV,1)); % and initialize it
        end % if
        out.(head{cfld})(idx_adcp)=dat{cfld+7}(idx_visea); % store visea data in output structure
    end
end


