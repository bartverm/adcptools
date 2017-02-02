classdef NMEA < handle
    properties
        raw_data
    end
    methods
        function obj=NMEA(varargin)
            narginchk(0,1);
            if nargin > 0
                obj.raw_data=varargin{1};
            end
        end
        function names=get_message_names(obj)
            names=fieldnames(obj.raw_data);                                               % Get the names of the available messages
        end
        function [msg_fids, fids]=get_fids(obj)
            msgs=obj.get_message_names();
            fids=[];
            msg_fids=cell(numel(msgs,1));
            for cf=1:numel(msgs)
                msg_fids{cf}=unique(obj.raw_data.(msgs{1}).file_id);
                fids=union(fids,msg_fids{cf});
            end
        end
        function [time, string_ids, file_ids]=get_time(obj)
            msgs=obj.get_message_names();                                               % Get the names of the available messages
            [msg_fids, fids]=obj.get_fids();
            time_field='utc';
            nfiles=numel(fids);
            [time, string_ids, file_ids]=deal(repmat({double.empty(0,1)}, nfiles,1));
            for cnt_msg=1:numel(msgs)
                c_msg=msgs{cnt_msg};
                if isfield(obj.raw_data.(c_msg),time_field)
                    for cnt_fid=1:numel(msg_fids{cnt_msg})
                        c_fid=msg_fids{cnt_msg}(cnt_fid);
                        tmptime=double(obj.raw_data.(c_msg).(time_field)(:,obj.raw_data.(c_msg).file_id==c_fid)');
                        tmptime=datenum(double([repmat([2000, 01, 01],size(tmptime,1),1), tmptime]));
                        time{fids==c_fid}=[time{fids==c_fid}; tmptime];
                        string_ids{fids==c_fid}=[string_ids{fids==c_fid} obj.raw_data.(c_msg).string_id(obj.raw_data.(c_msg).file_id==c_fid)];
                        file_ids{fids==c_fid}=[file_ids{fids==c_fid} obj.raw_data.(c_msg).file_id(obj.raw_data.(c_msg).file_id==c_fid)];
                    end
                end
            end
            time=vertcat(time{:})';
            string_ids=horzcat(string_ids{:});
            file_ids=horzcat(file_ids{:});
            fids=unique(file_ids);           
            for cnt_fid=1:numel(fids)
                cfid=fids(cnt_fid);
                filt=find(file_ids==cfid);
                [string_ids(filt), srt_idx]=sort(string_ids(filt));
                time(filt)=time(filt(srt_idx));
            end
            time=time+cumsum([0 diff(time)<-22/24]);
        end
    end
end