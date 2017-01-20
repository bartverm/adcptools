function NMEA=readNMEA(infiles, varargin)
%readNMEA(filenames) Reads in data from nmea text files
%   
%   NMEA=readNMEA(filenames) Reads NMEA data from files given in the
%       charachter array or cell of strings filenames. The output is a
%       structure, in which each field is one message type found in the
%       files. Each message is again a structure containing the data from
%       that message. For each message a string_id is given, which keeps
%       track of the line number the data was obtained from and a file_id
%       that keeps track of the file the data came from.
%
%   NMEA=readNMEA(filenames, 'ParamName', ParamValue) allows to set the
%       following parameters:
%           'RDens' 
%           {true} | false
%           Sets whether to read the non-standard RDENS strings
%
%           'Concatenate'
%           {true} | false
%           Sets whether to concatenate the output. If set to false it will
%           output an Nx1 cell array, with N being the number of input
%           files. Each cell contains a structure as described above, but
%           omitting the file_id information in the message structure
%
%   readNMEA currently supports the following standard messages:
%       $__GGA    : Position information. 
%       $__GLL    : Position information.
%       $__GSA    : Dilution of Precision (DOP) and active satellite
%                   information.
%       $__DBT    : Depth below transducer.
%       $__DBS    : Depth below surface.
%       $__ZDA    : Date and time information.
%       $__RMC    : Recommended minimum specific GPS data.
%       $__HDT    : Heading.
%       Standard NMEA formats are specified in the defineNMEA function.
%       Additional standard NMEA formats can be added in the defineNMEA
%       function.
%
%   readNMEA currently supports the following proprietary messages:
%       $RDENS    : RDI proprietary message containing ensemble number and 
%                   pc-time.
%
%   Author: Bart Vermeulen
%   Last edit: 21-12-2009

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

%% Handle input
P=inputParser;                                                             % create an inputParser object
P.addRequired('infiles',@(x) iscellstr(x) || ischar(x));                   % first input is the filenames
% P.addParameter('RDens',true,...
%     @(x) validateattributes(x,{'logical'},{'scalar'}));                    % Parameter to set whether to read RDENS messages
P.addParameter('Concatenate',true,...
    @(x) validateattributes(x, {'logical'},{'scalar'}));                   % Parameter to set whether to concatenate output
P.parse(infiles, varargin{:});                                             % Parse the input

% read_rdens=P.Results.RDens;                                                % Get value of RDens parameter
cat_cells=P.Results.Concatenate;                                           % Get value of Concatenate parameter

if ischar(infiles)                                                         % if one file is give as char
    infiles=cellstr(infiles);                                              % transform to scalar cell of char
end

% check if input files exist
for cf=length(infiles)                                                     % Loop over all filenames
    assert(exist(infiles{cf},'file')==2,...
        ['Couldn''t find file: ', infiles{cf}]);                           % Make sure file exists
end


%% Parse files
nfiles=numel(infiles);                                                     % get number of files
tmpNMEA=cell(size(infiles));                                               % initialize an arrray of structures

for cntfile=1:nfiles                                                       % loop over files
    disp(['Reading data from file ' infiles{cntfile}]);                    % display file being read
    fid=fopen(infiles{cntfile},'r');                                       % Open file
    if (fid < 0)                                                           % Check whether file was opened succesfully
        error('readNMEA:bad_fid','unable to open file for reading');       % If not, generate an error
    end
    
    lines=textscan(fid,'%s');                                              % Get all lines in the file
    fclose(fid);                                                           % Close the file
    lines=lines{1};                                                        % Take lines out of 1x1 cell
    
    tmpNMEA{cntfile}=parseNMEA(lines,true);                                % Parse the lines with NMEA parser
    
%     if read_rdens                                                          % If requested
%         tmpNMEA{cntfile}.rdens=parseRDENS(lines);                          % Parse the lines with RDENS parser
%     end
end % for cntfile

%% Concatenate outputs if requested
if cat_cells                                                               % If outputs must be concatenated
    fnames=cellfun(@fields, tmpNMEA, 'UniformOutput', false);              % Get all the message names
    fnames=unique(vertcat(fnames{:}));                                     % Retain only unique message names
    for c_field=1:numel(fnames)                                            % Loop over unique message names
        c_fname=fnames{c_field};                                           % Get current message name
        f_field=find(cellfun(@(x) isfield(x, c_fname), tmpNMEA));          % Find files which contain the current message (used also to store the file_id)
        tmpStruct=cellfun(@(x) x.(c_fname), tmpNMEA(f_field));             % Create N element structure with current message as field, where N is the number of files with the current message
        dnames=[fields(tmpStruct); {'file_id'}];                           % Get the names of the data fields in the current message and add the 'file_id' data to keep track of file
        for c_str=1:numel(tmpStruct)                                       % Loop over valid files
            tmpStruct(c_str).file_id=ones(1,size(...
                tmpStruct(c_str).(dnames{1}),2))*f_field(c_str);           % Create the file_id to keep track of originating file
        end
        for c_dat=1:numel(dnames)                                          % For each data field in current message
            c_dname=dnames{c_dat};                                         % Get the data field name                                  
            NMEA.(c_fname).(c_dname)=[tmpStruct.(c_dname)];                % Concatenate the data from each file
        end
    end
else                                                                       % If no concatenation is requested
    NMEA=tmpNMEA;                                                          % output data as is
end

