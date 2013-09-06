function adcpnmea = readNMEAADCP(inadcp,nmeafilename)
% READADCPNMEA reads data from files of the NMEA format and outputs
%              structures for each NMEA message type and its data order for
%              each ensemble. If multiple values are available for the same
%              ensemble these values are averaged.
%      
%              adcpnmea=readNMEAADCP(adcp,nmeafiles) reads data from the
%              files specified in the character array or cell of strings
%              NMEAFILES. The function matches data in the NMEA files with
%              using the ensemble numbers. If multiple ensembles are
%              present with the same number (This happens if you stopped
%              pinging and started again in the same data set) it assumes
%              files are given in chronological order.
%              When files contain overlapping messages, data is stored from
%              the last given file
%
%              Author: Bart Vermeulen
%              Last Edit: 23-08-2010
%              "[~," replaced by "[Dummy," for us in matlab 2009a (Frans)

%    Copyright 2009,2010 Bart Vermeulen, Frans Buschman
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

% Assumption: Files are given in chronological order and one series at a time!

inp=inputParser;                                                           % Create an object of the InputParser class
inp.addRequired('inadcp',@isstruct);                                       % Add required input inadcp
inp.addRequired('nmeafilename',@(x) (iscellstr(x) | ischar(x)));           % Add the required variable 'mmeafilename' and check for right format
inp.parse(inadcp,nmeafilename);                                            % Parse input
if ischar(nmeafilename)
    nmeafilename=cellstr(nmeafilename);                                    % Change character array to cell
end
clear inp

adcpnmea=struct;                                                           % Initialize adcpnmea as a structure
innmea=readNMEA(nmeafilename);                                             % Read the data from the NMEA files

if ~isfield(innmea,'RDENS')                                                % If no RDENS information is found
    error('readNMEAADCP:noRDENS','Cannot find Ensemble information')                 % generate an error
end
nmeamsgs=fieldnames(innmea);                                               % Get the names of the available messages

%% Find time and date corresponding to lineids

% Store times and dates and corresponding line ids
disp('Searching for time and date...')
linesTime=cell(size(innmea));                                               
UTCtimes=cell(size(innmea));
lines=cell(size(innmea));
hastime=false(size(innmea));
for cntmsg=1:length(nmeamsgs)                                              % Loop for all messages
    curmsg=nmeamsgs{cntmsg};                                               % store name of current message
    for cntfile=1:length(innmea);                                          % Loop for all files
        if isempty(innmea(cntfile).(curmsg))                               % If current message is empty in this file
            continue                                                       % Continue to next file
        end
        lines{cntfile}=[lines{cntfile};innmea(cntfile).(curmsg).lineid];   % Append to lines{current file} the lineids of this message
        if isfield(innmea(cntfile).(curmsg),'UTCtime')                     % If this message contains UTC time
            hastime(cntfile)=true;                                         % set flag for this file indicating this file contains time
            linesTime{cntfile}=[linesTime{cntfile};innmea(cntfile).(curmsg).lineid]; % Store the ids of the lines containing time for each file
            UTCtimes{cntfile}=[UTCtimes{cntfile};innmea(cntfile).(curmsg).UTCtime]; % Store the time corresponding to the stored lineids
        end
    end
end
clear curmsg cntmsg cntfile
% sort time and lines according to lines 
for cntfile=1:length(innmea)
    [linesTime{cntfile},idx]=unique(linesTime{cntfile});
    UTCtimes{cntfile}=UTCtimes{cntfile}(idx,:);
    lines{cntfile}=unique(lines{cntfile});
end


if any(hastime)
    %% Find time and date for all lineids
    timeV=cell(size(innmea));
    fhastime=find(hastime==1);
    curdayoffset=0;
    for cfile=1:length(fhastime)                                                       % Loop for all the files
        cntfile=fhastime(cfile);
        startdate=datenum([2008 1 1 0 0 0]);                               % Make up a date
        dayoffset=cumsum([curdayoffset;abs(diff(UTCtimes{cntfile}(:,1)))>22]);        % from jumps in hours deduce day jumps and store the day offset
        curdayoffset=dayoffset(end);
        tmpdate=repmat(startdate,[size(UTCtimes{cntfile},1),1]);           % Replicate startdate to match the amount of UTC time information
        tmpdate=datevec(tmpdate+dayoffset);                                % Add the day offset
        intY=datenum(double([tmpdate(:,1:3),UTCtimes{cntfile}]));          % Making datenumber from the timevector
        fnan=isnan(intY) | isnan(linesTime{cntfile});                      % find the nans 
        intY(fnan)=[];                                                     % remove them
        linesTime{cntfile}(fnan)=[];                                       % remove them
        [Dummy,idx]=unique(intY,'first');                                              
        while numel(idx)~=numel(intY)
            idx2=intersect(1:numel(intY),idx);
            intY(idx2)=intY(idx2)+1e-6;
            [Dummy,idx]=unique(intY,'first');
        end
        timeV{cntfile}=datevec(interp1(linesTime{cntfile},intY,...
                               lines{cntfile},'linear'));                  % interpolating time for lines without a time stamp
        lineswtime=find(all(~isnan(timeV{cntfile}),2));                    % Removing lines with no time 
        
        %% Find correspondence between time of ADCP and UTCtime
        fadcpfile=find(inadcp.FileNumber==cntfile);
        fline=(innmea(cntfile).RDENS.lineid >= lines{cntfile}(lineswtime(1)) & ...
            innmea(cntfile).RDENS.lineid <= lines{cntfile}(lineswtime(end)));         % search for RDENS lines with a timestamp
        [commens,idxnmea,idxadcp]=intersect(innmea(cntfile).RDENS.ensnum(fline),...
                                            inadcp.ensnum(fadcpfile)); % find all adcp ensemble with a corresponding RDENS with timestamp
        if isempty(commens)
            warning('readNMEAADCP2:noCorrespondance',['No RDENS within utc time fields in file: ',nmeafilename{cntfile}])
            hastime(cntfile)=0;
            continue
        elseif numel(commens)==1
            dt=datenum(indacp.timeV(fadcpfile(idxadcp),:))-datenum(timeV(lines{cntfile}==innmea(cntfile).RDENS.lineid(idxnmea),:));
        else
            dtstart=datenum(double(inadcp.timeV(fadcpfile(idxadcp(1)),:)))-datenum(timeV{cntfile}(lines{cntfile}==innmea(cntfile).RDENS.lineid(idxnmea(1)),:));
            dtend=datenum(double(inadcp.timeV(fadcpfile(idxadcp(end)),:)))-datenum(timeV{cntfile}(lines{cntfile}==innmea(cntfile).RDENS.lineid(idxnmea(end)),:));
            if abs(dtstart-dtend)>(5/24/60/60)
                warning('readNMEAADCP2:highDT','Time difference between adcp clock and UTCtime is changing too much, something might be wrong') 
            end
            dt=nanmean([dtstart,dtend]);
        end
       timeV{cntfile}=datevec(datenum(timeV{cntfile})+dt);
    end

    %% Find correspondance between data and timeV
    nmeamsgs(strcmp(nmeamsgs,'rdens'))=[];   % Remove RDENS
    tmpdata=struct;
    fhastime=find(hastime==1);
    for cntmsg = 1:length(nmeamsgs)                                            % Loop for all the available messages
        % Populate structure to contain all nmea data and corresponding
        % time of adcp
        actmsg=nmeamsgs{cntmsg};
        tmpdata.(actmsg)=struct;
        msgdata=fieldnames(innmea(find(cellfun(@isstruct,{innmea(:).(actmsg)}),1,'first')).(actmsg));                                % Get the available data for each message
        msgdata(strcmp(msgdata,'lineid'))=[];% Remove lineid from list of data availabe in the message
        tmpdata.(actmsg).timeV=[];
        for cntdat=1:length(msgdata)                                           % Loop for all the available data in a message
            actdt=msgdata{cntdat};
            tmpdata.(actmsg).(actdt)=[];
        end
        
        for cfile=1:length(fhastime)                                                       % Loop for all the files
            cntfile=fhastime(cfile);
            if isempty(innmea(cntfile).(actmsg))
                continue
            end
            [Dummy,idxa, idxb]=intersect(lines{cntfile},innmea(cntfile).(actmsg).lineid);
            ffile=isempty(tmpdata.(actmsg).timeV);
            if ~ffile
                lasttime=datenum(tmpdata.(actmsg).timeV(find(any(~isnan(tmpdata.(actmsg).timeV),2),1,'last'),:));
                tmpdata.(actmsg).timeV=[tmpdata.(actmsg).timeV;datevec(lasttime+0.0001/24/3600)];
            end
            tmpdata.(actmsg).timeV=[tmpdata.(actmsg).timeV;timeV{cntfile}(idxa,:)]; %interpolate time for each line from UTCtime
            for cntdat=1:length(msgdata)                                           % Loop for all the available data in a message
                actdt=msgdata{cntdat};
                if isempty(innmea(cntfile).(actmsg).(actdt))
                    continue
                end
                if ~ffile
                    tmpdata.(actmsg).(actdt)=[tmpdata.(actmsg).(actdt);tmpdata.(actmsg).(actdt)(end,:)*nan];
                end
                tmpdata.(actmsg).(actdt)=[tmpdata.(actmsg).(actdt);innmea(cntfile).(actmsg).(actdt)(idxb,:)];
            end
        end
    end

    %% Find data corresponding to ensembles
    for cntmsg = 1:length(nmeamsgs)
        actmsg=nmeamsgs{cntmsg};
%         outdata.(actmsg)=struct;
        msgdata=fieldnames(tmpdata.(actmsg));                                  % Get the available data for each message
        msgdata(strcmp(msgdata,'timeV'))=[];
        [xint,idx]=unique(datenum(tmpdata.(actmsg).timeV));
        nnan=~isnan(xint);
        xint=xint(nnan);
        idx=idx(nnan);

        for cntdat = 1:length(msgdata)
            actdt=msgdata{cntdat};
            if isempty(tmpdata.(actmsg).(actdt))
                continue
            end
            warning('off','MATLAB:interp1:NaNinY');
            if any(strcmp(actdt,{'heading','pitch','roll'}))
                adcpnmea.(actmsg).(actdt)=atan2(...
                    interp1(xint,sind(double(tmpdata.(actmsg).(actdt)(idx,:))),datenum(inadcp.timeV),'linear'),...
                    interp1(xint,cosd(double(tmpdata.(actmsg).(actdt)(idx,:))),datenum(inadcp.timeV),'linear'))/pi*180;            
            else
                adcpnmea.(actmsg).(actdt)=interp1(xint,double(tmpdata.(actmsg).(actdt)(idx,:)),datenum(inadcp.timeV),'linear');
            end
            warning('on','MATLAB:interp1:NaNinY');
        end
    end
else                                                           % If there is no time information use the old method
    warning('readNMEAADCP2:noTimeFound','Could not find time, using old method...')
    nmeamsgs(cat(1,cellfun(@(x) ~isempty(x),regexpi(nmeamsgs,'rdens'))))=[];   % Remove RDENS
    fnothastime=find(hastime==0);
    for cfiles = 1: length(fnothastime)                                 % Loop for all nmea files
        cntfiles=fnothastime(cfiles);
        for cntmsg = 1:length(nmeamsgs)                                            % Loop for all the available messages
            actmsg=nmeamsgs{cntmsg};
            if isempty(innmea(cntfiles).(actmsg))                          % If the message is empty in this file
                continue                                                   % Continue to next file
            end
            msgdata=fieldnames(innmea(cntfiles).(actmsg));                         % Get the available data for each message
            msgdata(cat(1,cellfun(@(x) ~isempty(x),regexpi(msgdata,'lineid'))))=[];% Remove lineid from list of data availabe in the message
            for cntdat=1:length(msgdata)                                           % Loop for all the available data in a message
                actdt=msgdata{cntdat};
                ensids=inadcp.ensnum(inadcp.FileNumber==cntfiles);                                                % Sort all available ensemble numbers
                for ensid = ensids                                             % Loop for each ensemble number in current adcp file
                    matchEnsADCP=inadcp.ensnum==ensid;                   % Find Ensemble numbers in ADCP file that match the given ensemble number
                    rdenstmp=[0;innmea(cntfiles).RDENS.lineid];                     % Add a zero to the rd ensemble line id's
                    frdens=find(innmea(cntfiles).RDENS.ensnum==ensid);
                    if isempty(frdens)                                          % If no RDENS string is available for the given ensemble number
                       continue                                                % Continue to the next ensemble
                    end
                    IsInEns=innmea(cntfiles).(actmsg).lineid>rdenstmp(frdens) & innmea(cntfiles).(actmsg).lineid<rdenstmp(frdens+1);               % Search which values are generated between this ensemble and the previous one
                    if any(strcmp(actdt,{'heading','pitch','roll'}))
                        adcpnmea.(actmsg).(actdt)(matchEnsADCP,:)=atan2(...
                        sind(nanmean(innmea(cntfiles).(actmsg).(actdt)(IsInEns,:),1)),...
                        cosd(nanmean(innmea(cntfiles).(actmsg).(actdt)(IsInEns,:),1)))/pi*180;% average all the data for this ensemble    
                    else
                        adcpnmea.(actmsg).(actdt)(matchEnsADCP,:)=...
                        nanmean(innmea(cntfiles).(actmsg).(actdt)(IsInEns,:),1);% average all the data for this ensemble    
                    end
                end    
            end
        end
    end
end
