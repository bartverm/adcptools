function dat = readViseaLogFiles(inadcp, fname)
% Reads VISEA log files and matches data with ADCP time
%
%   dat = read_visea_log_files(inadcp, fname) reads data in VISEA log
%   files. Currently only supports NMEA data that is read with the NMEA
%   toolbox. inadcp is the ADCP structure read with readADCP to which the
%   data are interpolated. fname defines the name of the log files to read.
%
%   see also: readDeployment, nmea

%    Copyright 2020 Bart Vermeulen
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

if exist('nmea.Message','class')~=8
    warning('Nmea toolbox not available, not reading VISEA log files')
    dat=struct;
    return;
end

visea_time_stamp_fmt='\[(?<year>\d{4})/(?<month>\d{2})/(?<day>\d{2}),\s(?<hour>\d{2}):(?<minute>\d{2}):(?<second>\d{2}\.\d{3})\]';

files = dir(fname);
folders={files(:).folder};
bytes=[files(:).bytes];
files={files(:).name};

% figure out different kinds of files
files_meta=regexp(files,'(?<dep_name>\w*)(?<number>\d{3})(?<name>(GPS|com\d{1,2})).log','names');
files_meta=[files_meta{:}];
all_types={files_meta(:).name};
f_types=unique(all_types);
timeadcp=datetime(inadcp.timeV);
for c_type=1:numel(f_types)
    cur_type=f_types{c_type};
    cur_idx=strcmp(all_types,cur_type);
    cur_files=files(cur_idx);
    cur_folders=folders(cur_idx);
    cur_bytes=bytes(cur_idx);
    rawdat(1,sum(cur_bytes))=' '; %#ok<AGROW>
    cumsiz=[0 cumsum(cur_bytes)];
    for cf=1:numel(cur_files)
        fid=fopen([cur_folders{cf},filesep, cur_files{cf}]);
        rawdat(cumsiz(cf)+1:cumsiz(cf+1))=fread(fid,'*char');
        fclose(fid);
    end

    % width of visea timestamp is always 26 characters
    [visea_t, indat, visea_pos, indat_pos]=regexp(rawdat,visea_time_stamp_fmt,'match','split','start','end');
    clearvars rawdat
    visea_t=[visea_t{:}];
    visea_t=textscan(visea_t,'[%f/%f/%f, %f:%f:%f]');
    visea_t=datetime(visea_t{:});

    % check if stamps are properly ordered?
    if seconds(diff(visea_t))<-1
        error('Files are not ordered in time')
    end

    % check what the remaining data look like



    indat_start=cumsum(cellfun(@numel,indat)); % this computes where each of the data chunks starts after concatenating chunks
    indat=[indat{:}]; % concatenate chunks of data
    indat_pos=indat_pos+1; % correct starting position of data chuncks in original file

    dat.(cur_type)=nmea.Message.parse_all(indat); % parse nmea strings
    if ~isempty(dat.(cur_type))
        fi=fields(dat.(cur_type));
        n_indat_start=numel(indat_start); % number of data chuncks
        for cf=1:numel(fi)
            cfn=fi{cf}; % current field name
            srt_raw=[reshape(indat_start,[],1);dat.(cur_type).(cfn).char_pos]; % concatenate position of chuncks and position of processed data 
            [~,i_raw]=sort(srt_raw); % order original chuncks position and data position
            f_dat_pos=find(i_raw>n_indat_start); % find where data is positioned among chunks
            pos_prev=1; % assume chunk index is one before
            pos_indat=max(i_raw(f_dat_pos-pos_prev),1); % get chunck number for each data field (this might end up being < 1!!!)
            in_excess=pos_indat > n_indat_start; % find data fields for which the assumption above was incorrect
            while any(in_excess)
                pos_prev=pos_prev+1; % assume chunk index is one more before
                pos_indat(in_excess)=max(i_raw(f_dat_pos(in_excess)-pos_prev),1); % get chunk number
                in_excess=pos_indat > n_indat_start; % check assumption again
            end
            pos_log=indat_pos(pos_indat); % find position of chunks in original file
            visea_time=reshape(interp1(visea_pos,visea_t,pos_log,'nearest','extrap'),[],1); % perform interpolation of vise time based on position in file
            % interpolate data to adcp time
            allf=fields(dat.(cur_type).(cfn));
            near_idx=interp1(visea_time,1:numel(visea_time),timeadcp,'nearest','extrap');
            for cfield=1:numel(allf)
                dat.(cur_type).(cfn).(allf{cfield})=dat.(cur_type).(cfn).(allf{cfield})(near_idx,:);
            end
        end
        clearvars rawdat
    end
    % parse pd0?
    
end
end
