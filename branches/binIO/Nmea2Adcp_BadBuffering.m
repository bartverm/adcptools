classdef Nmea2Adcp_BadBuffering < Nmea2Adcp
    % This class matches NMEA data to ADCP data in the case when the input data
    %   affected by poor handling of serial data flow in WinRiver.
    
    methods
        function obj=Nmea2Adcp_BadBuffering(varargin)
            obj@Nmea2Adcp(varargin{:});
        end
        function [string_id, file_id]=string_id_adcp(obj)
            adcp_ensnum=obj.pd0.get_ensnum();
            [string_id, file_id]=deal(nan(size(adcp_ensnum)));
            if ~isfield(nmea.raw_data,'rdens')
                error('No rdens message found')
            end
            
            [~, fids]=obj.get_fids();
            
            for cnt_file=1:numel(fids)
                cfid=fids(cnt_file);
                f_adcp_curfile=find(adcp_fileid==cfid);
                f_nmea_curfile=find(obj.nmea.raw_data.rdens.file_id==cfid);
                
                [~, pd0_idx, nmea_idx]=intersect(adcp_ensnum(f_adcp_curfile),obj.nmea.raw_data.rdens.ensnum(f_nmea_curfile));
                string_id(f_adcp_curfile(pd0_idx))=obj.nmea.raw_data.rdens.string_id(f_nmea_curfile(nmea_idx));
            end
            
        end
        function data=match(obj)
            
        end
        function idx=last_time(obj)
            max_lines_nmea=accumarray(nmea_file_ids',nmea_string_ids',[],@nanmax,nan,false);
            fgood_adcp=isfinite(adcp_file_id);
            max_lines_adcp=accumarray(adcp_file_id(fgood_adcp)',adcp_string_id(fgood_adcp)',[],@nanmax,nan,false);
            max_lines=max(max_lines_adcp,max_lines_nmea);
            
            [idx, time]=deal(nan(size(adcp_string_id)));
            for cnt_file=1:numel(max_lines)
                fline_ids=false(1,max_lines(cnt_file));
                fline_ids(nmea_string_ids(nmea_file_ids==cnt_file))=true;
                line_ids=double(fline_ids);
                line_ids(fline_ids)=find(fline_ids);
                % Could not think of a way to do the following in a vectorized way...damn it! Possibly slow
                last_good=0;
                for cl=1:numel(line_ids)
                    if line_ids(cl)==0
                        line_ids(cl)=last_good;
                    else
                        last_good=line_ids(cl);
                    end
                end
                % line_ids holds the last line in the file which contains time
                %     fline_ids=nmea_string_ids(line_ids(fgood));
                f_adcp_cur_file=find(adcp_file_id==cnt_file);
                f_nmea_cur_file=find(nmea_file_ids==cnt_file);
                % deal with case that rdens has no preceding zero
                lidx=line_ids(adcp_string_id(f_adcp_cur_file)); % Getting the line in the file of previous time for each adcp ensemble
                
                [~,nmea_idx_prev,adcp_idx_valid]=intersect(nmea_string_ids(f_nmea_cur_file),lidx);
                time(f_adcp_cur_file(adcp_idx_valid))=nmea_time(f_nmea_cur_file(nmea_idx_prev));
                idx(f_adcp_cur_file(adcp_idx_valid))=f_nmea_cur_file(nmea_idx_prev);
            end
        end
    end
end