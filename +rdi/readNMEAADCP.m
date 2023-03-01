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
inp.addRequired('nmeafilename',@(x) (isstring(x) | iscellstr(x) | ischar(x)));           % Add the required variable 'mmeafilename' and check for right format
inp.parse(inadcp,nmeafilename);                                            % Parse input
if ischar(nmeafilename)
    nmeafilename=cellstr(nmeafilename);                                    % Change character array to cell
end
clear inp

adcpnmea=struct;                                                           % Initialize adcpnmea as a structure
innmea=rdi.readNMEA(nmeafilename);                                             % Read the data from the NMEA files

if ~isfield(innmea,'RDENS')                                                % If no RDENS information is found
    error('readNMEAADCP:noRDENS','Cannot find Ensemble information')                 % generate an error
end
nmeamsgs=fieldnames(innmea);                                               % Get the names of the available messages

%% Find time and date corresponding to lineids

% Store times and dates and corresponding line ids
posTime=cell(size(innmea));                                               
UTCtimes=cell(size(innmea));
char_pos=cell(size(innmea));
hastime=false(size(innmea));
for cntmsg=1:length(nmeamsgs)                                              % Loop for all messages
    curmsg=nmeamsgs{cntmsg};                                               % store name of current message
    for cntfile=1:length(innmea)                                          % Loop for all files
        if isempty(innmea(cntfile).(curmsg))                               % If current message is empty in this file
            continue                                                       % Continue to next file
        end
        char_pos{cntfile}=[char_pos{cntfile};innmea(cntfile).(curmsg).char_pos];   % Append to lines{current file} the lineids of this message
        if isfield(innmea(cntfile).(curmsg),'utc')                         % If this message contains UTC time
            hastime(cntfile)=true;                                         % set flag for this file indicating this file contains time
            posTime{cntfile}=[posTime{cntfile};innmea(cntfile).(curmsg).char_pos]; % Store the ids of the lines containing time for each file
            UTCtimes{cntfile}=[UTCtimes{cntfile};innmea(cntfile).(curmsg).utc]; % Store the time corresponding to the stored lineids
        end
    end
end
clear curmsg cntmsg cntfile
% sort time and lines according to lines 
for cntfile=1:length(innmea)
    [posTime{cntfile},idx]=unique(posTime{cntfile});
    UTCtimes{cntfile}=UTCtimes{cntfile}(idx,:);
    char_pos{cntfile}=unique(char_pos{cntfile});
end


if any(hastime)
    %% Find time and date for all lineids
    pos_time=cell(size(innmea)); % time for each position in file
    fhastime=find(hastime==1);
    cur_day_offset=days(0);
    for cfile=1:length(fhastime)                                                       % Loop for all the files
        cntfile=fhastime(cfile);
        utc_dur = duration(UTCtimes{cntfile}); % get utc time as duration
        pos = posTime{cntfile};
        dayoffset = cumsum([cur_day_offset;...
            days(diff(utc_dur) < hours(-22))]);
        cur_day_offset=dayoffset(end);
        utc_dur = utc_dur + dayoffset;
        fnan=isnan(utc_dur) | isnan(pos);                      % find the nans 
        utc_dur(fnan)=[];                                                     % remove them
        pos(fnan)=[];                                       % remove them
        dt = diff(utc_dur);
        assert(~any(dt<0),'time is going backward in nmea data');

        % interpolate times to have all unique times
        f_zero_jump = [false; dt == 0];
        utc_dur(f_zero_jump) = interp1(...
            pos(~f_zero_jump),...
            utc_dur(~f_zero_jump),...
            pos(f_zero_jump));
        pos_time{cntfile}=interp1(pos, utc_dur,...
                               char_pos{cntfile},'linear');                  % interpolating time for lines without a time stamp
        lineswtime=find(~isnan(pos_time{cntfile}));                    % Removing lines with no time 
        
        %% Find correspondence between time of ADCP and UTCtime
        fadcpfile=find(inadcp.FileNumber==cntfile);
        if isempty(innmea(cntfile).RDENS)
            hastime(cntfile)=0;
            continue
        end
            
        fline=(innmea(cntfile).RDENS.char_pos >= char_pos{cntfile}(lineswtime(1)) & ...
            innmea(cntfile).RDENS.char_pos <= char_pos{cntfile}(lineswtime(end)));         % search for RDENS lines with a timestamp
        [commens,idxnmea,idxadcp]=intersect(innmea(cntfile).RDENS.ensnum(fline),...
                                            double(inadcp.ensnum(fadcpfile))); % find all adcp ensemble with a corresponding RDENS with timestamp
        if isempty(commens)
            warning('readNMEAADCP2:noCorrespondance',['No RDENS within utc time fields in file: ',nmeafilename{cntfile}])
            hastime(cntfile)=0;
            continue
        elseif numel(commens)==1
            dt=datetime(inadcp.timeV(fadcpfile(idxadcp),:))-pos_time(char_pos{cntfile}==innmea(cntfile).RDENS.char_pos(idxnmea));
        else
            dtstart=datetime(double(inadcp.timeV(fadcpfile(idxadcp(1)),:)))-pos_time{cntfile}(char_pos{cntfile}==innmea(cntfile).RDENS.char_pos(idxnmea(1)));
            dtend=datetime(double(inadcp.timeV(fadcpfile(idxadcp(end)),:)))-pos_time{cntfile}(char_pos{cntfile}==innmea(cntfile).RDENS.char_pos(idxnmea(end)));
            if abs(dtstart-dtend) > seconds(5)
                warning('readNMEAADCP2:highDT','Time difference between adcp clock and UTCtime is changing too much, something might be wrong') 
            end
            dt=mean([dtstart,dtend],'omitnan');
        end
       pos_time{cntfile}=pos_time{cntfile}+dt;
    end

    %% Find correspondance between data and timeV
    nmeamsgs(strcmpi(nmeamsgs,'rdens'))=[];   % Remove RDENS
    tmpdata=struct;
    fhastime=find(hastime==1);
    for cntmsg = 1:length(nmeamsgs)                                            % Loop for all the available messages
        % Populate structure to contain all nmea data and corresponding
        % time of adcp
        actmsg=nmeamsgs{cntmsg};
        tmpdata.(actmsg)=struct;
        msgdata=fieldnames(innmea(find(cellfun(@isstruct,{innmea(:).(actmsg)}),1,'first')).(actmsg));                                % Get the available data for each message
        msgdata(strcmp(msgdata,'char_pos'))=[];% Remove lineid from list of data availabe in the message
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
            [~,idxa, idxb]=intersect(char_pos{cntfile},innmea(cntfile).(actmsg).char_pos);
%             ffile=isempty(tmpdata.(actmsg).timeV); % first file (for appending)
%             if ~ffile
%                 lasttime=tmpdata.(actmsg).timeV(end);
%                 tmpdata.(actmsg).timeV=[tmpdata.(actmsg).timeV;datevec(lasttime+0.0001/24/3600)];
%             end
            tmpdata.(actmsg).timeV=[tmpdata.(actmsg).timeV;pos_time{cntfile}(idxa)]; %interpolate time for each line from UTCtime
            for cntdat=1:length(msgdata)                                           % Loop for all the available data in a message
                actdt=msgdata{cntdat};
                if isempty(innmea(cntfile).(actmsg).(actdt))
                    continue
                end
%                 if ~ffile
%                     tmpdata.(actmsg).(actdt)=[tmpdata.(actmsg).(actdt);tmpdata.(actmsg).(actdt)(end,:)*nan];
%                 end
                fempty = idxb > size(innmea(cntfile).(actmsg).(actdt),1);
                if isa(innmea(cntfile).(actmsg).(actdt),'datetime')
                    dat_in = NaT(size(idxb,1),size(innmea(cntfile).(actmsg).(actdt),2));
                    dat_in.TimeZone = innmea(cntfile).(actmsg).(actdt).TimeZone;
                else
                    dat_in = nan(size(idxb,1),size(innmea(cntfile).(actmsg).(actdt),2));
                end
                dat_in(~fempty ,:) = innmea(cntfile).(actmsg).(actdt)(idxb(~fempty),:);
                tmpdata.(actmsg).(actdt)=[tmpdata.(actmsg).(actdt);dat_in];
            end
        end
    end

    %%%%%% UPDATED TILL HERE!!
    %% Find data corresponding to ensembles
    for cntmsg = 1:length(nmeamsgs)
        actmsg=nmeamsgs{cntmsg};
%         outdata.(actmsg)=struct;
        msgdata=fieldnames(tmpdata.(actmsg));                                  % Get the available data for each message
        msgdata(strcmp(msgdata,'timeV'))=[];
        [xint,idx]=unique(tmpdata.(actmsg).timeV);
        nnat=~isnat(xint);
        xint=xint(nnat);
        idx=idx(nnat);

        for cntdat = 1:length(msgdata)
            actdt=msgdata{cntdat};
            if isempty(tmpdata.(actmsg).(actdt))
                continue
            end
            warning('off','MATLAB:interp1:NaNinY');
            if any(strcmp(actdt,{'heading','pitch','roll'}))
                adcpnmea.(actmsg).(actdt)=atan2d(...
                    interp1(xint,sind(double(tmpdata.(actmsg).(actdt)(idx,:))),datetime(inadcp.timeV),'linear'),...
                    interp1(xint,cosd(double(tmpdata.(actmsg).(actdt)(idx,:))),datetime(inadcp.timeV),'linear'));            
            else
                adcpnmea.(actmsg).(actdt)=interp1(xint,tmpdata.(actmsg).(actdt)(idx,:),datetime(inadcp.timeV),'linear');
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
            msgdata(cat(1,cellfun(@(x) ~isempty(x),regexpi(msgdata,'char_pos'))))=[];% Remove lineid from list of data availabe in the message
            for cntdat=1:length(msgdata)                                           % Loop for all the available data in a message
                actdt=msgdata{cntdat};
                ensids=inadcp.ensnum(inadcp.FileNumber==cntfiles);                                                % Sort all available ensemble numbers
                for ensid = ensids                                             % Loop for each ensemble number in current adcp file
                    matchEnsADCP=inadcp.ensnum==ensid;                   % Find Ensemble numbers in ADCP file that match the given ensemble number
                    rdenstmp=[0;innmea(cntfiles).RDENS.char_pos];                     % Add a zero to the rd ensemble line id's
                    frdens=find(innmea(cntfiles).RDENS.ensnum==ensid);
                    if isempty(frdens)                                          % If no RDENS string is available for the given ensemble number
                       continue                                                % Continue to the next ensemble
                    end
                    IsInEns=innmea(cntfiles).(actmsg).char_pos>rdenstmp(frdens) & innmea(cntfiles).(actmsg).char_pos<rdenstmp(frdens+1);               % Search which values are generated between this ensemble and the previous one
                    if any(strcmp(actdt,{'heading','pitch','roll'}))
                        adcpnmea.(actmsg).(actdt)(matchEnsADCP,:)=atan2d(...
                        sind(mean(innmea(cntfiles).(actmsg).(actdt)(IsInEns,:),1,'omitnan')),...
                        cosd(mean(innmea(cntfiles).(actmsg).(actdt)(IsInEns,:),1,'omitnan')));% average all the data for this ensemble    
                    else
                        adcpnmea.(actmsg).(actdt)(matchEnsADCP,:)=...
                        mean(innmea(cntfiles).(actmsg).(actdt)(IsInEns,:),1,'omitnan');% average all the data for this ensemble    
                    end
                end    
            end
        end
    end
end
