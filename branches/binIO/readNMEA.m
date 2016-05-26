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
P.addRequired('infiles',@(x) iscellstr(x) || ischar(x));
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
nfiles=length(infiles);   %find amount if fukes
NMEA(1:nfiles)=struct();  %initialize an arrray of structures

defineNMEA;

for cntfile=1:nfiles
    disp(['Reading data from file ' infiles{cntfile}]);
    fid=fopen(infiles{cntfile},'r');   %Open file
    if (fid < 0)
	error('readNMEA','unable to open file for reading');
    else
    fseek(fid,0,1);                    %Go to the end
    filesize=ftell(fid);               %Get filesize
    fseek(fid,0,-1);                   %Go to the beginning
    rawdat=fread(fid,filesize,'*char')'; %Read entire file
    fclose(fid);                       % Close the file
   
    
    
    % Search for NMEA and RDI strings
    [rawdat tok]=regexp(rawdat,patterns.nmea,'match','tokens');
    tok=vertcat(tok{:});

    if (size(tok,2) < 2)
        disp('Empty file or invalid NMEA data, skipping it');
    else
%     %Search for non RDENS data
%     csumidx=find(~ (strcmpi('RD',tok(:,1)) & strcmpi('ENS',tok(:,2))));
%     
%     %Do checksum test on all the rest
%     disp('Checking checksum...')
%     isvalnmea=nmeachecksum(rawdat(csumidx));
%     
%     %Remove all that failed checksum test
%     if ~all(isvalnmea)
%         disp(['Found ',num2str(length(find(~isvalnmea))),' sentences with bad checksum'])
%         rawdat(csumidx(~isvalnmea))=[];
%         tok(csumidx(~isvalnmea),:)=[];
%     end
    
    %Search for RDENS data
    rdensidx=strcmpi('RD',tok(:,1)) & strcmpi('ENS',tok(:,2));
    if any(rdensidx), disp('RDENS proprietary data found, reading...' )
        [NMEA(cntfile).RDENS,discard]=readRDENS(rawdat(rdensidx));
        NMEA(cntfile).RDENS.lineid=find(rdensidx);
        if ~isempty(discard)
            disp(['Discarding ',num2str(length(discard)),' malformed RDENS string(s)'])
            NMEA(cntfile).RDENS.lineid(discard)=[];
        end
    end
    clear rdensidx
     
    
    %Search for GGA data
    ggaidx=strcmpi('GGA',tok(:,2));
    if any(ggaidx), disp('GGA data found, reading...' )
         [NMEA(cntfile).GGA,discard]=readGGA(rawdat(ggaidx));
         NMEA(cntfile).GGA.lineid=find(ggaidx);
         if ~isempty(discard)
            disp(['Discarding ',num2str(length(discard)),' malformed GGA string(s)'])
            NMEA(cntfile).GGA.lineid(discard)=[];
         end
    end
    clear ggaidx
    
    %Search for GBS data
    gbsidx=strcmpi('PS',tok(:,1)) & strcmpi('AT',tok(:,2)) & strcmpi('GBS',tok(:,3));
    if any(gbsidx), disp('GBS proprietary data found, reading...' )
         [NMEA(cntfile).GBS,discard]=readPSAT(rawdat(gbsidx));
         NMEA(cntfile).GBS.lineid=find(gbsidx);
         if ~isempty(discard)
            disp(['Discarding ',num2str(length(discard)),' malformed GBS string(s)'])
            NMEA(cntfile).GBS.lineid(discard)=[];
         end
    end
    clear gbsidx
 
    %Search for GLL data
    gllidx=strcmpi('GLL',tok(:,2));
    if any(gllidx), disp('GLL data found, reading...' )
        [NMEA(cntfile).GLL,discard]=readGLL(rawdat(gllidx));
        NMEA(cntfile).GLL.lineid=find(gllidx);
        if ~isempty(discard)
           disp(['Discarding ',num2str(length(discard)),' malformed GLL string(s)'])
           NMEA(cntfile).GLL.lineid(discard)=[];
        end
    end
    clear gllidx

    %Search for GSA data
    gsaidx=strcmpi('GSA',tok(:,2));
    if any(gsaidx), disp('GSA data found, reading...' )
        [NMEA(cntfile).GSA,discard]=readGSA(rawdat(gsaidx));
        NMEA(cntfile).GSA.lineid=find(gsaidx);
        if ~isempty(discard)
            disp(['Discarding ',num2str(length(discard)),' malformed GSA string(s)'])
            NMEA(cntfile).GSA.lineid(discard)=[];
        end
    end
    clear gsaidx

    %Search for DBT data
    dbtidx=strcmpi('DBT',tok(:,2));
    if any(dbtidx), disp('DBT data found, reading...' )
        [NMEA(cntfile).DBT,discard]=readDBT(rawdat(dbtidx));
        NMEA(cntfile).DBT.lineid=find(dbtidx);
        if ~isempty(discard)
            disp(['Discarding ',num2str(length(discard)),' malformed DBT string(s)'])
            NMEA(cntfile).DBT.lineid(discard)=[];
        end
    end
    clear dbtidx

    % Search for HPR data
    hpridx=strcmpi('PS',tok(:,1)) & strcmpi('AT',tok(:,2)) & strcmpi('HPR',tok(:,3));
    if any(hpridx), disp('HPR data found, reading...' )
        [NMEA(cntfile).HPR,discard]=readHPR(rawdat(hpridx));
        NMEA(cntfile).HPR.lineid=find(hpridx);
        if ~isempty(discard)
            disp(['Discarding ',num2str(length(discard)),' malformed HPR string(s)'])
            NMEA(cntfile).HPR.lineid(discard)=[];
        end
    end
    clear hpridx
 
 
    %Search for DPT data (readDPT not yet implemented!)
%     dptidx=find(~(cellfun(@isempty,regexpi(rawdat,'\$..DPT'))));
%     if ~isempty(dptidx), disp('DPT data found, reading...' )
%         NMEA.DPT(cntfile)=readDPT(rawdat(dptidx));
%         NMEA.DPT(cntfile).lineid=dptidx';
%     end
%     clear dptidx

    %Search for DBS data
    dbsidx=strcmpi('DBS',tok(:,2));
    if any(dbsidx), disp('DBS data found, reading...' )
        [NMEA(cntfile).DBS,discard]=readDBS(rawdat(dbsidx));
        NMEA(cntfile).DBS.lineid=find(dbsidx);
        if ~isempty(discard)
            disp(['Discarding ',num2str(length(discard)),' malformed DBS string(s)'])
            NMEA(cntfile).DBS.lineid(discard)=[];
        end
    end
    clear dbsidx

    %Search for ZDA data
    zdaidx=strcmpi('ZDA',tok(:,2));
    if any(zdaidx), disp('ZDA data found, reading...' )
        [NMEA(cntfile).ZDA,discard]=readZDA(rawdat(zdaidx));
        NMEA(cntfile).ZDA.lineid=find(zdaidx);
        if ~isempty(discard)
            disp(['Discarding ',num2str(length(discard)),' malformed ZMA string(s)'])
            NMEA(cntfile).ZDA.lineid(discard)=[];
        end       
    end
    clear zdaidx

    %Search for RMC data
   rmcidx=strcmpi('RMC',tok(:,2));
     if any(rmcidx), disp('RMC data found, reading...' )
        [NMEA(cntfile).RMC,discard]=readRMC(rawdat(rmcidx));
        NMEA(cntfile).RMC.lineid=find(rmcidx);
        if ~isempty(discard)
            disp(['Discarding ',num2str(length(discard)),' malformed RMC string(s)'])
            NMEA(cntfile).RMC.lineid(discard)=[];
        end
    end
    clear rmcidx

    %Search for HDT data
    hdtidx=strcmpi('HDT',tok(:,2));
    if any(hdtidx), disp('HDT data found, reading...' )
        [NMEA(cntfile).HDT,discard]=readHDT(rawdat(hdtidx));
        NMEA(cntfile).HDT.lineid=find(hdtidx);
        if ~isempty(discard)
            disp(['Discarding ',num2str(length(discard)),' malformed HDT string(s)'])
            NMEA(cntfile).HDT.lineid(discard)=[];
        end
    end
    clear hdtidx
    
        % search for VTG data 
    vtgidx=strcmpi('VTG',tok(:,2));
    if any(vtgidx), disp('VTG data found, reading...' )
        [NMEA(cntfile).VTG,discard]=readVTG(rawdat(vtgidx));
        NMEA(cntfile).VTG.lineid=find(vtgidx);
        if ~isempty(discard)
            disp(['Discarding ',num2str(length(discard)),' malformed VTG string(s)'])
            NMEA(cntfile).VTG.lineid(discard)=[];
        end
    end
    clear vtgidx
%     
% 
% %     gstpos=regexpi(rawdat,'\$..GST');
% %     rrepos=regexpi(rawdat,'\$..RRE');
% 
% %     hdgpos=regexpi(rawdat,'\$..HDG');
% %     hdmpos=regexpi(rawdat,'\$..HDM');
% %     rotpos=regexpi(rawdat,'\$..ROT');
    end % else of if size(tok,2) < 2
    end % else of if (fid < 0) 
end % for cntfile

