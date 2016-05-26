function dataout=readADCP(files,varargin)
%READADCP Reads a RDI adcp binary data file
%         [ADCP]=readADCP(FILENAME) reads adcp data and outputs all data
%         
%         FILENAME can be a character array with filenames on each row or a
%         cell of strings containing filenames in each element. The
%         function will read multiple files even when configuration differs
%         amongst them. In this case the function will issue a warning.
%         Settings will be stored once for each file
%
%         [ADCP]=readADCP(FILENAME,VARS) reads adcp data given in VARS
%         VARS should be a character array containing one of the following
%         characters:
%         e   ensemble info
%         v   velocity data
%         c   correlation data
%         h   echo intensity
%         p   percentage good
%         b   bottom track data
%         x   external data stored in the adcp file
%         f   external data stored in external winriver ascii-files (The
%             function will look for the files based on the raw-data
%             filenames. Being implemented!
%
%         Example: ADCP=readADCP('data_000.000; data_001.000','vEyx') will
%         read velocity, ensemble information and external data from the
%         two files
%
%         Author:      Bart Vermeulen, David Vermaas, Maximiliano Sassi
%         Last edit:   26-03-2009

%    Copyright 2007-2010 Bart Vermeulen, David Vermaas, Maximiliano Sassi
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


WRII_GGA_ID=104;
WRII_HDT_ID=107;
WRII_VTG_ID=101;
WRII_DBT_ID=102;

% global variables are: data_defs, data_out, ens, data, raw_dat


%% Check input
if verLessThan('matlab','7.13')
    error(nargchk(1,2,nargin)) %#ok<NCHKN>
    error(nargoutchk(0,1,nargout)) %#ok<NCHKE>
else
    narginchk(1,2)
    nargoutchk(0,1)
end
inp=inputParser;                                                           % Create an object of the InputParser class
inp.addRequired('files',@(x) (iscellstr(x) | ischar(x)));                  % Add the required variable 'files' and check for right format
inp.addOptional('flags','vhcpbex',@ischar);                                % Add the optional argument 'flags' and check for right format
inp.parse(files,varargin{:});                                              % Parse input

if ischar(files)
    files=cellstr(files);                                                % Change character array to cell
end
flags=inp.Results.flags;
clear inp

%Check filenames
nfiles=numel(files);
raw_dat=cell(nfiles,1);
filesize=zeros(nfiles,1);
disp('Reading files');
for cntfile=1:nfiles                                                       % Loop for files
    [fid,openmessage]=fopen(files{cntfile},'r');                           % Open file in read-only binary mode
    if fid==-1                                                             % If opening file fails
        warning('readADCP:WrongFile',[openmessage,': ',files{cntfile},...
            ', Skipping this file.'])                                      % Show warning
        continue                                                           % Continue to next file
    end
    [raw_dat{cntfile}, filesize(cntfile)]=fread(fid,'*uint8');             % read all data in file
    fclose(fid);                                                           % close file
end
raw_dat=cat(1,raw_dat{:});                                                 % concatenate raw data

% search for ensembles
ens.pos=find_valid_ensembles();
assert(~isempty(ens.pos),'readADCP:NoEnsemble','Could not find any valid ensemble')      % Generate error

% locate the data fields
ens.ndat=double(raw_dat(ens.pos+5)');
data=locate_data_fields();

disp(['Found ',num2str(numel(ens.pos)), ' ensembles and ', num2str(numel(data.ensid)),' data fields.']);

% parse the fixed leader
read_data_fields()



%% Preallocate and read data



% % Read external data
% if any(regexpi(flags,'x'))                                                 % If external data is to be read
%     if any(AllHeaders(:,1)==8226)                                          % Check if any WinRiver II data is available
%         if any(AllHeaders(:,2)==WRII_GGA_ID)                               % If WRII GGA data is available
%             disp('WinRiver II GGA data found')                             % Display that this data is found
%             dd=diff(find(AllHeaders(:,2)==WRII_GGA_ID &...
%                 AllHeaders(:,1)==8226));                                   % Determine distance between all WRII GGA headers
%             nblocks=max(diff(find(dd~=1)));                                % Determine maximum distance (in datablocks) between non WRII GGA headers (i.e. maximum WRII GGA headers in one ensembles)
%             initNMEAGGA(nens,nblocks);                                     % Initialize variables for WRII GGA data
%             for cntens=1:nens                                              % Loop for every ensemble
%                 Ndatablock=find(DataHeader{cntens}(:,2)==WRII_GGA_ID &...
%                     DataHeader{cntens}(:,1)==8226);                        % Determine Number(s) of datablock(s) with WRII GGA data
%                 if isempty(Ndatablock)                                     % If no WRII GGA data is available
%                     continue                                               % Continue with next ensemble
%                 else                                                       % Otherwise
%                     for cntDatBlock=1:length(Ndatablock)                   % Loop for every available datablock
%                         fpos=EnsStart{cntens}+...
%                             DataOffset{cntens}(Ndatablock(cntDatBlock));   % Determine position of data
%                         readNMEAGGA(fileid(cntens),fpos,...
%                             cntens,cntDatBlock);                           % Call function to read out the data
%                     end
%                 end
%             end
%         end
%         if any(AllHeaders(:,2)==WRII_HDT_ID)                               % Same as above but for WRII HDT data
%             disp('WinRiver II HDT data found')
%             dd=diff(find(AllHeaders(:,2)==WRII_HDT_ID...
%                 & AllHeaders(:,1)==8226));
%             nblocks=max(diff(find(dd~=1)));
%             initNMEAHDT(nens,nblocks);
%             for cntens=1:nens
%                 Ndatablock=find(DataHeader{cntens}(:,2)==WRII_HDT_ID &...
%                     DataHeader{cntens}(:,1)==8226);
%                 if isempty(Ndatablock)
%                     continue
%                 else
%                     for cntDatBlock=1:length(Ndatablock)
%                         fpos=EnsStart{cntens}+...
%                             DataOffset{cntens}(Ndatablock(cntDatBlock));
%                         readNMEAHDT(fileid(cntens),fpos,cntens,cntDatBlock);
%                     end
%                 end
%             end
%         end
%         if any(AllHeaders(:,2)==WRII_VTG_ID)                               % Same as above but for WRII VTG data
%             disp('WinRiver II VTG data found')
%             dd=diff(find(AllHeaders(:,2)==WRII_VTG_ID &...
%                 AllHeaders(:,1)==8226));
%             nblocks=max(diff(find(dd~=1)));
%             initNMEAVTG(nens,nblocks);
%             for cntens=1:nens
%                 Ndatablock=find(DataHeader{cntens}(:,2)==WRII_VTG_ID &...
%                     DataHeader{cntens}(:,1)==8226);
%                 if isempty(Ndatablock)
%                     continue
%                 else
%                     for cntDatBlock=1:length(Ndatablock)
%                         fpos=EnsStart{cntens}+...
%                             DataOffset{cntens}(Ndatablock(cntDatBlock));
%                         readNMEAVTG(fileid(cntens),fpos,cntens,cntDatBlock);
%                     end
%                 end
%             end
%         end
%         if any(AllHeaders(:,2)==WRII_DBT_ID)                               % Same as above but for WRII DBT data
%             disp('WinRiver II DBT data found')
%             dd=diff(find(AllHeaders(:,2)==WRII_DBT_ID &...
%                 AllHeaders(:,1)==8226));
%             nblocks=max(diff(find(dd~=1)));
%             initNMEADBT(nens,nblocks);
%             for cntens=1:nens
%                 Ndatablock=find(DataHeader{cntens}(:,2)==WRII_DBT_ID &...
%                     DataHeader{cntens}(:,1)==8226);
%                 if isempty(Ndatablock)
%                     continue
%                 else
%                     for cntDatBlock=1:length(Ndatablock)
%                         fpos=EnsStart{cntens}+...
%                             DataOffset{cntens}(Ndatablock(cntDatBlock));
%                         readNMEADBT(fileid(cntens),fpos,cntens,cntDatBlock);
%                     end
%                 end
%             end
%         end
%     end
%     if any(AllHeaders(:,1)==8448)                                           % Same as Velocity data but for WinRiver DBT data
%         tmpdbt=cell(nens,1);
%         for cntens=1:nens
%             Ndatablock=find(DataHeader{cntens}(:,1)==8448);
%             if isempty(Ndatablock)
%                 continue
%             else
%                 fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
%                 tmpdbt{cntens}=getDBT(fileid(cntens),fpos);
%             end
%         end
%         readDBTint(tmpdbt,nens);
%         clear tmpdbt;    
%     end
%     if any(AllHeaders(:,1)==8449)                                          % Sama as above but for WR GGA data
%         tmpgga=cell(nens,1);
%         for cntens=1:nens
%             Ndatablock=find(DataHeader{cntens}(:,1)==8449);
%             if isempty(Ndatablock)
%                 continue
%             else
%                 fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
%                 tmpgga{cntens}=getGGA(fileid(cntens),fpos);
%             end
%         end
%         readGGAint(tmpgga,nens);
%         clear tmpgga
%     end
%     if any(AllHeaders(:,1)==8450)                                          % Same as above but for WR VTG data
%         tmpvtg=cell(nens,1);
%         for cntens=1:nens
%             Ndatablock=find(DataHeader{cntens}(:,1)==8450);
%             if isempty(Ndatablock)
%                 continue
%             else
%                 fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
%                 tmpvtg{cntens}=getVTG(fileid(cntens),fpos);
%             end
%         end
%         readVTGint(tmpvtg,nens);
%         clear tmpvtg;
%     end
%     if any(AllHeaders(:,1)==8452)                                          % Same as above but for WR HDT data
%         tmphdt=cell(nens,1);
%         for cntens=1:nens
%             Ndatablock=find(DataHeader{cntens}(:,1)==8452);
%             if isempty(Ndatablock)
%                 continue
%             else
%                 fpos=EnsStart{cntens}+DataOffset{cntens}(Ndatablock);
%                 tmphdt{cntens}=getHDT(fileid(cntens),fpos);
%             end
%         end
%         readHDTint(tmphdt,nens);
%         clear tmphdt;
%     end
% end



%% READING FUNCTIONS

function read_data_fields()
    data_defs=define_data_types();
    for count_data=1:numel(data_defs)
        switch data_defs(count_data).type
            case 'fields'
                for count_field=1:numel(data_defs(count_data).fields)
                    current_field=data_defs(count_data).fields(count_field);
                    dataout.(current_field.name)=get_scalar_data(data_defs(count_data).header,current_field);
                end
            case 'array'
                if ~isfield(dataout,'nbins') || ~isfield(dataout,'usedbeams')
                    warning('readADCP:NoVariableLeader','Cannot read array data(velocity, echo, etc) without variable leader data')
                    continue
                end
                current_field=data_defs(count_data).fields(1);
                dataout.(current_field.name)=get_array_data(data_defs(count_data).header,current_field);
            case 'nmea'
        end
    end
end



%Read velocity
function V=readVEL(fid,fpos,nbins)
    fseek(fid,fpos,-1);
    dataID=fread(fid,2,'*uint8');
    if any(dataID~=[0;1])
        error('readADCP:readVEL:wrongID','Velocity data ID seems to be wrong')
    end
    V=reshape(fread(fid,4*double(nbins),'*int16'),4,[])';
end

%Read echo intensity
function H=readECHO(fid,fpos,nbins)
    fseek(fid,fpos,-1);
    dataID=fread(fid,2,'*uint8');
    if any(dataID~=[0;3])
        error('readADCP:readECHO:wrongID','Echo intesity data ID seems to be wrong')
    end
    H=reshape(fread(fid,4*double(nbins),'*uint8'),4,[])';    
end

%Read correlation magnitude
function C=readCORR(fid,fpos,nbins)
    fseek(fid,fpos,-1);
    dataID=fread(fid,2,'*uint8');
    if any(dataID~=[0;2])
        error('readADCP:readCORR:wrongID','Correlation data ID seems to be wrong')
    end
    C=reshape(fread(fid,4*double(nbins),'*uint8'),4,[])';    
end

%Read percentage good
function P=readPERC(fid,fpos,nbins)
    fseek(fid,fpos,-1);
    dataID=fread(fid,2,'*uint8');
    if any(dataID~=[0;4])
        error('readADCP:readPERC:wrongID','Percentage good data ID seems to be wrong')
    end
    P=reshape(fread(fid,4*double(nbins),'*uint8'),4,[])';    
end

% Read bottom track data
function readBT(fid,fpos,cntens)
    fseek(fid,fpos,-1);
    dataID=fread(fid,2,'*uint8');
    if any(dataID~=[0;6])
        error('readADCP:readBT:wrongID','Bottom-Track data ID seems to be wrong')
    end
    bt1=fread(fid,2,'*uint16');
    dataout.btpingperens(cntens)=bt1(1);                  %Bottom tracking pings per ensemble (0-999)
    dataout.reacqdelay(cntens)=bt1(2);                    %Delay in number of ensembles before reacquiring (0-999)
    bt2=fread(fid,4,'*uint8');
    dataout.mincormag(cntens)=bt2(1);                      %Minimum for correlation magnitude in counts (0-255)
    dataout.minevampl(cntens)=bt2(2);                      %Minimum evaluation amplitude in counts (1-255)
    dataout.btminpergood(cntens)=bt2(3);                   %Minimum percentage of good bt pings 
    dataout.btmode(cntens)=bt2(4);                         %BT mode
    dataout.btmaxerrv(cntens)=fread(fid,1,'*uint16');                     %Maximum bt error velocity in mm/s (0-5000)
    fseek(fid,4,0);                                                 %Reserved data
    bt3=fread(fid,4,'uint16=>uint32');
    dataout.btrange(cntens,1:4)=bt3(1:4);                      %bt range of beam 1,2,3,4 in cm (0-65535)
    bt4=fread(fid,4,'*int16');
    dataout.btvel(cntens,1:4)=bt4(1:4);                         %bt velocty of beam 1,2,3,4 in mm/s (-32768 to 32768)
    bt5=fread(fid,12,'*uint8');
    dataout.btcor(cntens,1:4)=bt5(1:4);                         %bt correlation magnitude beam 1,2,3,4 in counts (0-255)
    dataout.btevampl(cntens,1:4)=bt5(5:8);                      %bt evaluation amplitude for strength of bottom echo beam 1,2,3,4 in counts (0-255)
    dataout.btpercgood(cntens,1:4)=bt5(9:12);                    %bt percentage of good pings in beam 1,2,3,4 (0-100)
    bt6=fread(fid,3,'*uint16');
    dataout.minlyrsize(cntens)=bt6(1);                    %Minimum size of ref layer in dm (0-9999)
    dataout.nearbnd(cntens)=bt6(2);                       %Near boundary of ref layer in dm (0-9999)
    dataout.farbnd(cntens)=bt6(3);                        %Far boundary of ref layer in dm (0-9999) 
    bt7=fread(fid,4,'*int16');    
    dataout.reflyrvel(cntens,1:4)=bt7(1:4);                     %Reference layer velocity 1,2,3,4 in mm/s (-32768 to 32768)
    bt8=fread(fid,12,'*uint8');
    dataout.reflyrcor(cntens,1:4)=bt8(1:4);                     %Reference layer correlation magnitude beam 1,2,3,4 in counts (0-255)
    dataout.reflyrint(cntens,1:4)=bt8(5:8);                     %Reference layer echo intensity beam 1,2,3,4 in counts (0-255)
    dataout.reflyrpergood(cntens,1:4)=bt8(9:12);                 %Reference layer percentage good pings in beam 1,2,3,4 (0-100)
    dataout.maxdepth(cntens)=fread(fid,1,'*uint16');                      %bt maximum depth in dm (80-9999)
    bt9=fread(fid,5,'*uint8');
    dataout.rssiamp(cntens,1:4)=bt9(1:4);                       %Received signal strength indicator in beam 1,2,3,4 in counts (0-255), 1 count ca. 0.45 dB
    dataout.gain(cntens)=bt9(5);                           %Gain level for shallow water
    bt10=fread(fid,4,'uint8=>uint32');
    dataout.btrange(cntens,1:4)=reshape(bt10(1:4)*65636,[1,4])+dataout.btrange(cntens,1:4);     %Most significant byte of bt range in beam 1,2,3,4 in cm (65636-16777215)
end


% Read WinRiverII General NMEA GGA data
function readNMEAGGA(fid,fpos,cntens,cntblock)
    fseek(fid,fpos,-1);
    dataID=fread(fid,2,'*uint16');
    if any(dataID~=[8226;WRII_GGA_ID])
        error('readADCP:readNMEAGGA:wrongID','WinRiver II general NMEA GGA data ID seems to be wrong')
    end
    fseek(fid,2,0);
    dataout.NMEAGGA.deltaT(cntens,cntblock)=fread(fid,1,'*float64');
    dataout.NMEAGGA.msgHeader{cntens}(cntblock,:)=fread(fid,7,'uchar=>char')';
    dataout.NMEAGGA.UTC{cntens}(cntblock,:)=fread(fid,10,'uchar=>char')';
    dataout.NMEAGGA.Lat(cntens,cntblock)=fread(fid,1,'*float64');
    dataout.NMEAGGA.SN(cntens,cntblock)=fread(fid,1,'uchar=>char');
    dataout.NMEAGGA.Long(cntens,cntblock)=fread(fid,1,'*float64');
    dataout.NMEAGGA.EW(cntens,cntblock)=fread(fid,1,'uchar=>char');
    dataout.NMEAGGA.Qual(cntens,cntblock)=fread(fid,1,'*uint8');
    dataout.NMEAGGA.NSat(cntens,cntblock)=fread(fid,1,'*uint8');
    dataout.NMEAGGA.HDOP(cntens,cntblock)=fread(fid,1,'*float32');
    dataout.NMEAGGA.Alt(cntens,cntblock)=fread(fid,1,'*float32');
    dataout.NMEAGGA.AltUnit(cntens,cntblock)=fread(fid,1,'uchar=>char');
    dataout.NMEAGGA.Geoid(cntens,cntblock)=fread(fid,1,'*float32');
    dataout.NMEAGGA.GeoidUnit(cntens,cntblock)=fread(fid,1,'uchar=>char');
    dataout.NMEAGGA.AgeDGPS(cntens,cntblock)=fread(fid,1,'*float32');
    dataout.NMEAGGA.RefStatID(cntens,cntblock)=fread(fid,1,'*uint16');
end

% Read WinRiverII General NMEA HDT data
function readNMEAHDT(fid,fpos,cntens,cntblock)
    fseek(fid,fpos,-1);
    dataID=fread(fid,2,'*uint16');
    if any(dataID~=[8226;WRII_HDT_ID])
        error('readADCP:readNMEAHDT:wrongID','WinRiver II general NMEA HDT data ID seems to be wrong')
    end
    fseek(fid,2,0);
    dataout.NMEAHDT.deltaT(cntens,cntblock)=fread(fid,1,'*float64');
    dataout.NMEAHDT.msgHeader{cntens}(cntblock,:)=fread(fid,7,'uchar=>char')';
    dataout.NMEAHDT.heading(cntens,cntblock)=fread(fid,1,'*float64');
    dataout.NMEAHDT.trueInd(cntens,cntblock)=fread(fid,1,'uchar=>char');
end

% Read WinRiverII General NMEA VTG data
function readNMEAVTG(fid,fpos,cntens,cntblock)
    fseek(fid,fpos,-1);
    dataID=fread(fid,2,'*uint16');
    if any(dataID~=[8226;WRII_VTG_ID])
        error('readADCP:readNMEAVTG:wrongID','WinRiver II general NMEA VTG data ID seems to be wrong')
    end
    fseek(fid,2,0);
    dataout.NMEAVTG.deltaT(cntens,cntblock)=fread(fid,1,'*float64');
    dataout.NMEAVTG.msgHeader{cntens}(cntblock,:)=fread(fid,7,'uchar=>char')';
    dataout.NMEAVTG.COGTrue(cntens,cntblock)=fread(fid,1,'*float64');
    dataout.NMEAVTG.TrueIndicator(cntens,cntblock)=fread(fid,1,'uchar=>char');
    dataout.NMEAVTG.COGMagn(cntens,cntblock)=fread(fid,1,'*float64');
    dataout.NMEAVTG.MagnIndicator(cntens,cntblock)=fread(fid,1,'uchar=>char');
    dataout.NMEAVTG.SpeedOverGroundKts(cntens,cntblock)=fread(fid,1,'*float64');
    dataout.NMEAVTG.KtsIndicator(cntens,cntblock)=fread(fid,1,'uchar=>char');
    dataout.NMEAVTG.SpeedOverGroundKmh(cntens,cntblock)=fread(fid,1,'*float64');
    dataout.NMEAVTG.KmhIndicator(cntens,cntblock)=fread(fid,1,'uchar=>char');
    dataout.NMEAVTG.ModeIndicator(cntens,cntblock)=fread(fid,1,'uchar=>char');
end

% Read WinRiverII General NMEA DBT data
function readNMEADBT(fid,fpos,cntens,cntblock)
    fseek(fid,fpos,-1);
    dataID=fread(fid,2,'*uint16');
    if any(dataID~=[8226;WRII_DBT_ID])
        error('readADCP:readNMEADBT:wrongID','WinRiver II general NMEA DBT data ID seems to be wrong')
    end
    fseek(fid,2,0);
    dataout.NMEADBT.deltaT(cntens,cntblock)=fread(fid,1,'*float64');
    dataout.NMEADBT.msgHeader{cntens}(cntblock,:)=fread(fid,7,'uchar=>char')';
    dataout.NMEADBT.WaterDepthFt(cntens,cntblock)=fread(fid,1,'*float64');
    dataout.NMEADBT.FeetIndicator(cntens,cntblock)=fread(fid,1,'uchar=>char');
    dataout.NMEADBT.WaterDepthm(cntens,cntblock)=fread(fid,1,'*float64');
    dataout.NMEADBT.MeterIndicator(cntens,cntblock)=fread(fid,1,'uchar=>char');
    dataout.NMEADBT.WaterDepthF(cntens,cntblock)=fread(fid,1,'*float64');
    dataout.NMEADBT.FathomIndicator(cntens,cntblock)=fread(fid,1,'uchar=>char');
end

% Read WinRiver GGA data
function ggastr=getGGA(fid,fpos)
    fseek(fid,fpos,-1);
    dataID=fread(fid,2,'*uint8');
    if any(dataID~=[1;33])
        error('readADCP:getGGA:wrongID','WinRiver GGA data ID seems to be wrong')
    end
    ggastr=fgetl(fid);
end

function readGGAint(ggastr,nens)
    defineNMEA;
    hasdata=~cellfun(@isempty,regexp(ggastr,patterns.gga));
    if ~any(hasdata)
        return
    end
    disp('WinRiver GGA data found')
    initGGA(nens);
    datgga=readGGA(ggastr(hasdata));
    dataout.GGA.UTCtime(hasdata,:)=datgga.UTCtime;
    dataout.GGA.lat(hasdata)=datgga.lat;
    dataout.GGA.long(hasdata)=datgga.long;
    dataout.GGA.qualind(hasdata)=datgga.qualind;
    dataout.GGA.numsat(hasdata)=datgga.numsat;
    dataout.GGA.hdop(hasdata)=datgga.hdop;
    dataout.GGA.antalt(hasdata)=datgga.antalt;
    dataout.GGA.geosep(hasdata)=datgga.geosep;
    dataout.GGA.agediff(hasdata)=datgga.agediff;
    dataout.GGA.diffid(hasdata)=datgga.diffid;
end

%Read WinRiver HDT data
function hdtstr=getHDT(fid,fpos)
    fseek(fid,fpos,-1);
    dataID=fread(fid,2,'*uint8');
    if any(dataID~=[4;33])
        error('readADCP:getHDT:wrongID','WinRiver HDT data ID seems to be wrong')
    end
    hdtstr=fgetl(fid);
end

function readHDTint(hdtstr,nens)
    defineNMEA;
    hasdata=~cellfun(@isempty,regexp(hdtstr,patterns.hdt));
    if ~any(hasdata)
        return
    end
    disp('WinRiver HDT data found')
    initHDT(nens);
    dathdt=readHDT(hdtstr(hasdata));
    dataout.HDT.heading(hasdata)=dathdt.heading;
end

%Read WinRiver DBT data
function dbtstr=getDBT(fid,fpos)
    fseek(fid,fpos,-1);
    dataID=fread(fid,2,'*uint8');
    if any(dataID~=[0;33])
        error('readADCP:getDBT:wrongID','WinRiver DBT data ID seems to be wrong')
    end
    dbtstr=fgetl(fid);
end

function readDBTint(dbtstr,nens)
    defineNMEA;
    hasdata=~cellfun(@isempty,regexp(dbtstr,patterns.dbt));
    if ~any(hasdata)
        return
    end
    disp('WinRiver DBT data found')
    initDBT(nens);
    datdbt=readDBT(dbtstr(hasdata));
    dataout.DBT.depthf(hasdata)=datdbt.depthf;
    dataout.DBT.dephtM(hasdata)=datdbt.depthM;
    dataout.DBT.depthF(hasdata)=datdbt.depthF;
end

%Read WinRiver VTG data
function vtgstr=getVTG(fid,fpos)
    fseek(fid,fpos,-1);
    dataID=fread(fid,2,'*uint8');
    if any(dataID~=[2;33])
        error('readADCP:getVTG:wrongID','WinRiver VTG data ID seems to be wrong')
    end
    vtgstr=fgetl(fid);
end

function readVTGint(vtgstr,nens)
    defineNMEA
    hasdata=~cellfun(@isempty,regexp(vtgstr,patterns.vtg));
    if ~any(hasdata)
        return
    end
    disp('WinRiver VTG data found')
    initVTG(nens);
    datvtg=readVTG(vtgstr(hasdata));
    dataout.VTG.TrackDegTrue(hasdata)=datvtg.TrackDegTrue;
    dataout.VTG.TrackDegMagn(hasdata)=datvtg.TrackDegMagn;
    dataout.VTG.SpeedKnots(hasdata)=datvtg.SpeedKnots;
    dataout.VTG.SpeedKmH(hasdata)=datvtg.SpeedKmH;
    dataout.VTG.mode(hasdata)=datvtg.mode;
end

%%          PREALLOCATION FUNCTIONS

%Initialize Variable leader
function initVL(nens)
    dataout.ensnum= zeros(1,nens, 'uint32');
    dataout.BITcheck= char(zeros(nens,8));
    dataout.speedsound= zeros(1,nens, 'uint16');
    dataout.depthtransd= zeros(1,nens, 'uint16');
    dataout.heading= zeros(1,nens, 'uint16');
    dataout.pitch= zeros(1,nens, 'int16');
    dataout.roll= zeros(1,nens, 'int16');
    dataout.salinity= zeros(1,nens, 'uint16');
    dataout.temperature= zeros(1,nens, 'int16');
    dataout.prepingT= zeros(1,nens, 'uint16');
    dataout.headstd= zeros(1,nens, 'uint8');
    dataout.pitchstd= zeros(1,nens, 'uint8');
    dataout.rollstd= zeros(1,nens, 'uint8');
    dataout.ADC= zeros(nens,8, 'uint8');
    dataout.errorstat1= zeros(nens,8, 'uint8');
    dataout.errorstat2= zeros(nens,8, 'uint8');
    dataout.errorstat3= zeros(nens,8, 'uint8');
    dataout.errorstat4= zeros(nens,8, 'uint8');
    dataout.pressure= zeros(1,nens, 'uint32');
    dataout.pressurevar= zeros(1,nens, 'uint32');
    dataout.timeV= zeros(nens,6, 'double');
end

%Initialize Bottom tracking data
function initBT(nens)
    dataout.btpingperens=zeros(1,nens,'uint16');
    dataout.reacqdelay=zeros(1,nens,'uint16');
    dataout.mincormag=zeros(1,nens,'uint8');
    dataout.minevampl=zeros(1,nens,'uint8');
    dataout.btminpergood=zeros(1,nens,'uint8');
    dataout.btmode=zeros(1,nens,'uint8');
    dataout.btmaxerrv=zeros(1,nens,'uint16');
    dataout.btrange=zeros([nens,4],'uint32');
    dataout.btvel=zeros([nens,4],'int16');
    dataout.btcor=zeros([nens,4],'uint8');
    dataout.btevampl=zeros([nens,4],'uint8');
    dataout.btpercgood=zeros([nens,4],'uint8');
    dataout.minlyrsize=zeros(1,nens,'uint16');
    dataout.nearbnd=zeros(1,nens,'uint16');
    dataout.farbnd=zeros(1,nens,'uint16');
    dataout.reflyrvel=zeros([nens,4],'int16');
    dataout.reflyrcor=zeros([nens,4],'uint8');
    dataout.reflyrint=zeros([nens,4],'uint8');
    dataout.reflyrpergood=zeros([nens,4],'uint8');
    dataout.maxdepth=zeros(1,nens,'uint16');
    dataout.rssiamp=zeros([nens,4],'uint8');
    dataout.gain=zeros(1,nens,'uint8');
end

%Initialize WinRiverII external GGA data
function initNMEAGGA(nens,nblocks)
    dataout.NMEAGGA.deltaT=ones(nens,nblocks,'double')*NaN;
    dataout.NMEAGGA.msgHeader=cell(nens,1);      %char(ones(nens,7,'uint8'));
    dataout.NMEAGGA.UTC=cell(nens,1);            %char(ones(nens,10,'uint8'));
    dataout.NMEAGGA.Lat=ones(nens,nblocks,'double')*NaN;
    dataout.NMEAGGA.SN=char(ones(nens,nblocks,'uint8'));
    dataout.NMEAGGA.Long=ones(nens,nblocks,'double')*NaN;
    dataout.NMEAGGA.EW=char(ones(nens,nblocks,'uint8'));
    dataout.NMEAGGA.Qual=ones(nens,nblocks,'uint8');
    dataout.NMEAGGA.NSat=ones(nens,nblocks,'uint8');
    dataout.NMEAGGA.HDOP=ones(nens,nblocks,'single');
    dataout.NMEAGGA.Alt=ones(nens,nblocks,'single')*NaN;
    dataout.NMEAGGA.AltUnit=char(ones(nens,nblocks,'uint8'));
    dataout.NMEAGGA.Geoid=ones(nens,nblocks,'single')*NaN;
    dataout.NMEAGGA.GeoidUnit=char(ones(nens,nblocks,'uint8'));
    dataout.NMEAGGA.AgeDGPS=ones(nens,nblocks,'single')*NaN;
    dataout.NMEAGGA.RefStatID=ones(nens,nblocks,'uint16');
end

%Initialize WinRiverII external GGA data
function initNMEAHDT(nens,nblocks)
    dataout.NMEAHDT.deltaT=ones(nens,nblocks,'double')*NaN;
    dataout.NMEAHDT.msgHeader=cell(nens,1);
    dataout.NMEAHDT.heading=ones(nens,nblocks,'double')*NaN;
    dataout.NMEAHDT.trueInd=char(ones(nens,nblocks,'uint8'));
end

%Initialize WinRiverII external GGA data
function initNMEAVTG(nens,nblocks)
    dataout.NMEAVTG.deltaT=ones(nens,nblocks,'double')*NaN;
    dataout.NMEAVTG.msgHeader=cell(nens,1);
    dataout.NMEAVTG.COGTrue=ones(nens,nblocks,'double')*NaN;
    dataout.NMEAVTG.TrueIndicator=char(ones(nens,nblocks,'uint8'));
    dataout.NMEAVTG.COGMagn=ones(nens,nblocks,'double')*NaN;
    dataout.NMEAVTG.MagnIndicator=char(ones(nens,nblocks,'uint8'));
    dataout.NMEAVTG.SpeedOverGroundKts=ones(nens,nblocks,'double')*NaN;
    dataout.NMEAVTG.KtsIndicator=char(ones(nens,nblocks,'uint8'));
    dataout.NMEAVTG.SpeedOverGroundKmh=ones(nens,nblocks,'double')*NaN;
    dataout.NMEAVTG.KmhIndicator=char(ones(nens,nblocks,'uint8'));
    dataout.NMEAVTG.ModeIndicator=char(ones(nens,nblocks,'uint8'));
end

%Initialize WinRiverII external GGA data
function initNMEADBT(nens,nblocks)
    dataout.NMEADBT.deltaT=ones(nens,nblocks,'double')*NaN;
    dataout.NMEADBT.msgHeader=cell(nens,1);
    dataout.NMEADBT.WaterDepthFt=ones(nens,nblocks,'double')*NaN;
    dataout.NMEADBT.FeetIndicator=char(ones(nens,nblocks,'uint8'));
    dataout.NMEADBT.WaterDepthm=ones(nens,nblocks,'double')*NaN;
    dataout.NMEADBT.MeterIndicator=char(ones(nens,nblocks,'uint8'));
    dataout.NMEADBT.WaterDepthF=ones(nens,nblocks,'double')*NaN;
    dataout.NMEADBT.FathomIndicator=char(ones(nens,nblocks,'uint8'));
end

%Initialize WinRiver DBT external data
function initDBT(nens)
    dataout.DBT.depthf=nan(1,nens,'single');
    dataout.DBT.dephtM=nan(1,nens,'single');
    dataout.DBT.depthF=nan(1,nens,'single');
end

%Initialize WinRiver GGA external data
function initGGA(nens)
    dataout.GGA.UTCtime=nan(nens,3,'single');
    dataout.GGA.lat=nan(1,nens,'double');
    dataout.GGA.long=nan(1,nens,'double');
    dataout.GGA.qualind=zeros(1,nens,'uint8');
    dataout.GGA.numsat=zeros(1,nens,'uint8');
    dataout.GGA.hdop=nan(1,nens,'single');
    dataout.GGA.antalt=nan(1,nens,'single');
    dataout.GGA.geosep=nan(1,nens,'single');
    dataout.GGA.agediff=nan(1,nens,'single');
    dataout.GGA.diffid=zeros(1,nens,'uint16');
end

%Initialize WinRiver VTG external data
function initVTG(nens)
    dataout.VTG.TrackDegTrue=nan(1,nens,'single');
    dataout.VTG.TrackDegMagn=nan(1,nens,'single');
    dataout.VTG.SpeedKnots=nan(1,nens,'single');
    dataout.VTG.SpeedKmH=nan(1,nens,'single');
    dataout.VTG.mode=zeros(1,nens,'uint8');
end

%Initialize WinRiver HDT external data
function initHDT(nens)
    dataout.HDT.heading=nan(1,nens,'single');
end

    %%          FUNCTION TO CHECK DIFFERENCES IN FIXED LEADER
    function issame=CheckFL(nValidFiles)
    issame=true;                                                               % Initialize result as true
    if nValidFiles<2                                                           % If less than two files
        issame=true;                                                           % Return true
        return                                                                 % Exit function (Nothing to check!!)
    end
    for cntfiles3=2:nValidFiles                                                 % Loop for all files
        issame=issame &&...                                                    % Check if any of the variables changes
            isequal(dataout.firmver(1),dataout.firmver(cntfiles3))&&...
        isequal(dataout.firmrev(1),dataout.firmrev(cntfiles3))&&...
        isequal(dataout.sysconf(1,:),dataout.sysconf(cntfiles3,:))&&...
        isequal(dataout.SymData(1),dataout.SymData(cntfiles3))&&...
        isequal(dataout.LagLength(1),dataout.LagLength(cntfiles3))&&...
        isequal(dataout.usedbeams(1),dataout.usedbeams(cntfiles3))&&...
        isequal(dataout.nbins(1),dataout.nbins(cntfiles3))&&...
        isequal(dataout.pingperens(1),dataout.pingperens(cntfiles3))&&...
        isequal(dataout.binsize(1),dataout.binsize(cntfiles3))&&...
        isequal(dataout.blnk(1),dataout.blnk(cntfiles3))&&...
        isequal(dataout.minthrsh(1),dataout.minthrsh(cntfiles3))&&...
        isequal(dataout.ncodrep(1),dataout.ncodrep(cntfiles3))&&...
        isequal(dataout.minpercgood(1),dataout.minpercgood(cntfiles3))&&...
        isequal(dataout.maxerrvel(1),dataout.maxerrvel(cntfiles3))&&...
        isequal(dataout.Tbetweenpng(1),dataout.Tbetweenpng(cntfiles3))&&...
        isequal(dataout.corinfo(1,:),dataout.corinfo(cntfiles3,:))&&...
        isequal(dataout.headalign(1),dataout.headalign(cntfiles3))&&...
        isequal(dataout.headbias(1),dataout.headbias(cntfiles3))&&...
        isequal(dataout.sensource(1,:),dataout.sensource(cntfiles3,:))&&...
        isequal(dataout.senavail(1,:),dataout.senavail(cntfiles3,:))&&...
        isequal(dataout.distmidbin1(1),dataout.distmidbin1(cntfiles3))&&...
        isequal(dataout.lngthtranspulse(1),dataout.lngthtranspulse(cntfiles3))&&...
        isequal(dataout.watrefbins(1,:),dataout.watrefbins(cntfiles3,:))&&...
        isequal(dataout.mintarget(1),dataout.mintarget(cntfiles3))&&...
        isequal(dataout.lowlattrig(1),dataout.lowlattrig(cntfiles3))&&...
        isequal(dataout.distpulse(1),dataout.distpulse(cntfiles3))&&...
        isequal(dataout.cpuserial(1,:),dataout.cpuserial(cntfiles3,:))&&...
        isequal(dataout.bandwidth(1),dataout.bandwidth(cntfiles3))&&...
        isequal(dataout.syspower(1),dataout.syspower(cntfiles3))&&...
        isequal(dataout.basefreqid(1),dataout.basefreqid(cntfiles3))&&...
        isequal(dataout.serial(1,:),dataout.serial(cntfiles3,:))&&...
        isequal(dataout.HADCPbeamangle(1),dataout.HADCPbeamangle(cntfiles3));
        if issame==false
            break
        end
    end
    end

    function data_idx=find_last_datafield(header)
        data_idx=nan(1,numel(ens.pos));
        ftype=find(data.headers==header);
        data_idx(data.ensid(ftype))=ftype; % Keeps last, if multiple data blocks of same type. nan means data not found in ensemble
    end

    function val=get_scalar_data(header,field_meta)
        dataidx=find_last_datafield(header);
        fgood=isfinite(dataidx);
        if ~any(fgood), val=[]; return, end;
        nval=nullval(field_meta.type);
        val=nval(ones(field_meta.size,numel(ens.pos)));
        val(:,data.ensid(dataidx(fgood)))=parse_blocks(bsxfun(@plus,(0:field_meta.size-1)'* sizeof(field_meta.type),data.offset(dataidx(fgood))+field_meta.offset-1),field_meta.type);
    end

    function val=get_array_data(header, field_meta)
        nval=nullval(field_meta.type);
        nbins_max=double(max(dataout.nbins));
        nbeams_max=double(max(dataout.usedbeams));
        val=nval(ones(nbins_max,numel(ens.pos),nbeams_max)); % Initialize
        if isempty(val), return, end;

        ndat=double(dataout.nbins).*double(dataout.usedbeams);
        % Index computations below only need to be performed once, now they are done on every array read action!
        
        ensidx=zeros(1,sum(ndat)); 
        ensidx(cumsum(ndat(1:end-1))+1)=1;
        ensidx=cumsum(ensidx)+1;

        % Beam indices
        cnvels=cumsum([0  ndat]);
        beamidx=mod(cumsum(ones(1,sum(ndat)))-cnvels(ensidx)-1,double(dataout.usedbeams(ensidx)))+1; % This is slow

        % Cell indices
        cncells=cumsum([0 double(dataout.nbins)]);
        cellidx=cumsum(beamidx==1)-cncells(ensidx);

        idces=sub2ind([double(max(dataout.nbins)),numel(ens.pos),double(max(dataout.usedbeams))],cellidx,ensidx,beamidx);

        % Velocity position
        pos=find_last_datafield(header); % Locate all  data
            % account for possible ensembles without data (does this ever happen?)
            f_bad_ens=find(~isfinite(pos));
            [~,f_bad_idx,~]=intersect(ensidx,f_bad_ens); % Get index of all velocities which are missing
            ensidx(f_bad_idx)=[];
            ndat(f_bad_ens)=[];
        idces(f_bad_idx)=[];
        pos=data.offset(pos(ensidx))+2+cumsum([0 sizeof(field_meta.type)*ones(1,sum(ndat)-1)]);
        cnvels=cumsum([0 ndat]);
        pos=pos-cnvels(ensidx)*sizeof(field_meta.type);

        val(idces)=parse_blocks(pos,field_meta.type);
    end

    function data=locate_data_fields()
        n_data=sum(ens.ndat);
        data.ensid=zeros(1,n_data);
        data.ensid(cumsum(ens.ndat)+1)=1;
        data.ensid=cumsum(data.ensid)+1;
        data.ensid(end)=[];
        pos_data_offset=(1:n_data);
        tmp_n_data_in_ens=[0 cumsum(ens.ndat)];
        pos_data_offset=pos_data_offset-tmp_n_data_in_ens(data.ensid);
        data.offset=double(parse_blocks(ens.pos(data.ensid)+6+(pos_data_offset-1)*2,'uint16'))+ens.pos(data.ensid);
        data.headers=parse_blocks(data.offset,'uint16');
    end

    function ens_pos=find_valid_ensembles()
        ens_pos=[];
        if numel(raw_dat)<4, return, end

        % Search for valid ensemble headers
        ens_pos=find([all([raw_dat(1:end-3)==127 raw_dat(2:end-2)==127],2); false])'; % find valid ensemble headers
        ens_bytes=parse_blocks(ens_pos+2,'uint16');
        if isempty(ens_pos), return, end

        % Remove if header position + ensbytes exceeds size of buffer
        fbad=ens_pos+double(ens_bytes)-1>numel(raw_dat)-2; 
        ens_pos(fbad)=[]; 
        ens_bytes(fbad)=[]; 
        if isempty(ens_pos), return, end


        % Remove if checksum in file does not match with computed one
        csum=[0 cumsum(double(raw_dat))']; 
        csum=uint16(mod(csum(ens_pos+double(ens_bytes))-csum(ens_pos),65536));
        csumr=typecast(reshape(raw_dat(bsxfun(@plus,ens_pos+double(ens_bytes),[0;1])),[],1),'uint16')';
        fbad=csumr~=csum;
        ens_pos(fbad)=[]; 
        ens_bytes(fbad)=[]; 
        if isempty(ens_pos), return, end

        % Remove if ens_pos + ens_bytes +2 is higher than start of next
        % ensemble (Ensembles shouldn't overlap, but don't need to be
        % contiguous), and yes, this happens even after checksum checks.
        while numel(ens_pos)>1 % loop while there's still more than 2 ensembles
            fbad=find((ens_pos(1:end-1)+double(ens_bytes(1:end-1))+2)>ens_pos(2:end))+1; % Search for overlapping ensembles
            if isempty(fbad), break, end % if no bad ensemble is found, exit loop
            fbad=fbad([true; diff(fbad)~=1]); % remove only first of series of consecutive overlapping ensembles
            ens_pos(fbad)=[]; % remove
            ens_bytes(fbad)=[]; % remove
        end
    end

    function out=parse_blocks(pos,type)
        if isempty(pos), out = []; return, end

        if strcmp(type,'uint8')
            out=raw_dat(pos);
        else
            nb=sizeof(type);
            out=typecast(reshape(raw_dat(bsxfun(@plus,pos,(0:nb-1)')),[],1),type)';
        end
    end
end


%%          INTERPRATATION FUNCTIONS

% Function to interpret system id
function sysstr=sysinterp(sysid)
    sysstrtemp='';
    switch sysid(1:3)
        case '000'
            sysstrtemp=[sysstrtemp,'75 KHz System'];
        case '100'
            sysstrtemp=[sysstrtemp,'150 KHz System'];
        case '010'
            sysstrtemp=[sysstrtemp,'300 KHz System'];
        case '110'
            sysstrtemp=[sysstrtemp,'600 KHz System'];
        case '001'
            sysstrtemp=[sysstrtemp,'1200 KHz System'];
        case '101'
            sysstrtemp=[sysstrtemp,'2400 KHz System'];
    end
    if sysid(4)=='0'
        sysstrtemp=[sysstrtemp,', Concave beam pattern'];
    else
        sysstrtemp=[sysstrtemp,', Convex beam pattern'];
    end
    switch sysid(5:6)
        case '00'
        sysstrtemp=[sysstrtemp,', Sensor configuration #1'];
        case '10'
        sysstrtemp=[sysstrtemp,', Sensor configuration #2'];
        case '11'
        sysstrtemp=[sysstrtemp,', Sensor configuration #3'];
    end
    if sysid(7)=='0'
        sysstrtemp=[sysstrtemp,', XDCR HD Not attached'];
    else
        sysstrtemp=[sysstrtemp,', XDCR HD Attached'];
    end
    if sysid(8)=='0'
        sysstrtemp=[sysstrtemp,', Up facing beam'];
    else
        sysstrtemp=[sysstrtemp,', Down facing beam'];
    end
    switch sysid(9:10)
        case '00'
        sysstrtemp=[sysstrtemp,', Beam angle: 15E'];
        case '10'
        sysstrtemp=[sysstrtemp,', Beam angle: 20E'];
        case '11'
        sysstrtemp=[sysstrtemp,', Beam angle: 30E'];
        case '01'
        sysstrtemp=[sysstrtemp,', Unknown beam angle'];
    end
    switch sysid(13:16)
        case '0010'
            sysstrtemp=[sysstrtemp,', 4 Beam Janus configuration'];
        case '1010'
            sysstrtemp=[sysstrtemp,', 5 Beam Janus configuration, 3 Demods'];
        case '1111'
            sysstrtemp=[sysstrtemp,', 5 Beam Janus configuration, 2 Demods'];
    end
sysstr=sysstrtemp;
end

% Function to interprete coordinate information (bits are reversed wrt manual of ADCP)
function corstr=corinterpr(corid)
            corstrtemp='';
            switch corid(4:5)
                case '00'
                    corstrtemp=[corstrtemp,'Beam Coordinates Used'];
                case '10'
                    corstrtemp=[corstrtemp,'Instrument Coordinates Used'];
                case '01'
                    corstrtemp=[corstrtemp,'Ship Coordinates Used'];
                case '11'
                    corstrtemp=[corstrtemp,'Earth Coordinates Used'];
            end
            if corid(1)=='1' 
                corstrtemp=[corstrtemp,', Bin mapping used'];
            else
                corstrtemp=[corstrtemp,', Bin mapping not used'];                  
            end
            if corid(2)=='1'
                corstrtemp=[corstrtemp,', 3 beams solution allowed'];
            else
                corstrtemp=[corstrtemp,', 3 beams solution not allowed'];                  
            end
            if corid(3)=='1'
                corstrtemp=[corstrtemp,', Tilts used (only for ship or eart coordinates)'];
            else
                corstrtemp=[corstrtemp,', Tilts not used'];                  
            end
corstr=corstrtemp;
end

function valens=isvalens(fid,headpos)
% Check validity of ensemble

fmes=fseek(fid,headpos+2,-1);                                              %Go to byte after two header bytes
if fmes==-1                                                                %Check if move was possible
    valens=false;
    return
end
EnsBytes=fread(fid,1,'uint16=>double');                                    %read number of bytes in ensemble
fmes=fseek(fid,headpos+EnsBytes+2,-1);                                     %Go to the end of the ensemble to check if file doesn't end before
if fmes==-1
    valens=false;
    return
end
fseek(fid,headpos,-1);                                                     %move to beginning of ensemble 
sumbytes=dec2bin(sum(fread(fid,EnsBytes,'uint8')));                        %calculate checksum
checksum=fread(fid,1,'*uint16');                                           %read checksum
if checksum==bin2dec(sumbytes(max((length(sumbytes)-15),1):end))                    %Check checksum
    valens=true;
else
    valens=false;
end
end

function headpos=searchHead(fid,fpos)
% search in the file for a valid header
headpos=-1;
fmes=fseek(fid,fpos,-1);
if fmes==-1
    return
end
while (~feof(fid))&& headpos==-1                                           %Loop until a end of file or until a header is found
    head=fread(fid,1,'*uint8');                                            %Read 1st byte
    if isempty(head) || head~=127                                          %Check for valid header value
        continue                                                           %If bad continue searching
    end
    head=fread(fid,1,'*uint8');                                            %Read 2nd byte
    if isempty(head)|| head~=127                                           %Check for valid header value
        continue                                                           %If bad continue searching
    end
    headpos=ftell(fid)-2;
end
end

% Function to check ensemble validity and read pos of ensemble and its data
function [NDataTypes,DataOffset,DataHeader]=readhead(fid,headpos)
fseek(fid,headpos+5,-1);                                                   %Move to the 'number of data-types'-field
NDataTypes=fread(fid,1,'uint8=>double');                                   %Read the number of data-types
DataOffset=fread(fid,NDataTypes,'uint16');                                 %Read offsets to data-types
DataHeader=ones(NDataTypes,2,'uint16')*65535;                              %Initialize header vector
for cntDataType=1:NDataTypes
    fseek(fid,headpos+DataOffset(cntDataType),-1);                         %Move to data header in ensemble
    DataHeader(cntDataType,:)=fread(fid,2,'*uint16')';                     %read header of data type (And next byte for navigation data)
end
end





function nb=sizeof(type)
    nb=cast(1,type); %#ok<NASGU>
    nb=whos('nb');
    nb=nb.bytes;
end





function val=nullval(type)
    if strcmp(type(1:3),'int'), val=intmin(type); 
    elseif strcmp(type(1:4),'uint'), val = intmax(type);
    elseif strcmp(type,{'single','double'}), val=nan;
    else val=0;
    end
end