function NMEA=readNMEA(infiles)
%readNMEA(filenames) Reads in data from nmea text files
%   
%   NMEAstruct=readNMEA(filenames)Reads NMEA data from files given in the
%       charachter array or cell of strings filenames and returns a 
%       Nfiles x 1 data structure containing one structure for each NMEA 
%       message encountered. For description of the several NMEA messages 
%       refer to the corresponding reading functions (i.e. for a $_GGA 
%       message refer to the help of readGGA.
%
%   readNMEA supports the following messages:
%       $RDENS    : RDI proprietary message containing ensemble number and 
%                   pc-time. Read by readRDENS
%       $PSAT,GBS : CSI proprietary message containing error estimates of
%       $PSAT,HPR : CSI proprietary message containing heading pitch and
%                   roll. Read by readHPR
%       $__GGA    : Position information. Read by readGGA
%                   GPS position. Read by readGBS
%       $__GLL    : Position information. Read by readGLL
%       $__GSA    : Dilution of Precision (DOP) and active satellite
%                   information. Read by readGSA
%       $__DBT    : Depth below transducer. Read by readDBT
%       $__DBS    : Depth below surface. Read by readDBS
%       $__ZDA    : Date and time information. Read by readZDA
%       $__RMC    : Recommended minimum specific GPS data. Read by readRMC
%       $__HDT    : Heading. Read by readHDT
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
P=inputParser;
P.addRequired('infiles',@(x) iscellstr(x) || isstring(x) || ischar(x));
P.parse(infiles);
infiles=P.Results.infiles;

if ischar(infiles)
    infiles=cellstr(infiles);
end

% check if input files exist
for cf=length(infiles)
    assert(exist(infiles{cf},'file')==2,['Couldn''t find file: ', infiles{cf}]);
end

%% Parse files
nfiles=length(infiles);   %find number of files
% NMEA(1:nfiles)=struct();  %initialize an arrray of structures

for cntfile=1:nfiles
    fid=fopen(infiles{cntfile},'r');   %Open file
    if (fid < 0)
	    error('readNMEA:CouldNotOpenFile','unable to open file for reading');
    end
    rawdat=fread(fid,'*char')'; %Read entire file
    fclose(fid);                       % Close the file
   
    NMEA(cntfile) = nmea.Message.parse_all(rawdat);
    
end % for cntfile

