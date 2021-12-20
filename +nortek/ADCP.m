classdef ADCP < ADCP
    properties(Constant, Access = protected)
        BURST_ID = 0x15
        AVERAGE_ID = 0x16
        BOTTOM_TRACKING_ID = 0x17
        INTERLEAVED_BURST_ID = 0x18
        BURST_ALTIMETER_ID = 0x1A
        DVL_BOTTOM_TRACKING_ID = 0x1B
        ECHO_SOUNDER_ID = 0x1C
        DVL_WATER_TRACKING_ID = 0x1D
        ALTIMETER_ID = 0x1E
        AVERAGE_ALTIMETER_ID = 0x1F
        STRING_ID = 0xA0
    end
    properties(Dependent, SetAccess = protected)
        nbytes
        has_burst
        has_average
        has_bottom_tracking
        has_interleaved_burst
        has_burst_altimeter
        has_dvl_bottom_tracking
        has_echo_sounder
        has_dvl_water_tracking
        has_altimeter
        has_average_altimeter
        has_string
    end
    methods
        function obj = ADCP(varargin)
            obj = obj@ADCP(varargin{:})
            for ca = 1 : nargin
                if ischar(varargin{ca}) || isStringScalar(varargin{ca})
                    obj.read_file(varargin{ca})
                end
            end
        end
        function val = get.nbytes(obj)
            val = numel(obj.raw);
        end
        function val = get.has_burst(obj)
            val = any(obj.data_type == obj.BURST_ID);
        end
        function val = get.has_average(obj)
            val = any(obj.data_type == obj.AVERAGE_ID);
        end
        function val = get.has_bottom_tracking(obj)
            val = any(obj.data_type == obj.BOTTOM_TRACKING_ID);
        end
        function val = get.has_interleaved_burst(obj)
            val = any(obj.data_type == obj.INTERLEAVED_BURST_ID);
        end
        function val = get.has_burst_altimeter(obj)
            val = any(obj.data_type == obj.BURST_ALTIMETER_ID);
        end
        function val = get.has_dvl_bottom_tracking(obj)
            val = any(obj.data_type == obj.DVL_BOTTOM_TRACKING_ID);
        end
        function val = get.has_echo_sounder(obj)
            val = any(obj.data_type == obj.ECHO_SOUNDER_ID);
        end
        function val = get.has_dvl_water_tracking(obj)
            val = any(obj.data_type == obj.DVL_WATER_TRACKING_ID);
        end
        function val = get.has_altimeter(obj)
            val = any(obj.data_type == obj.ALTIMETER_ID);
        end
        function val = get.has_average_altimeter(obj)
            val = any(obj.data_type == obj.AVERAGE_ALTIMETER_ID);
        end
        function val = get.has_string(obj)
            val = any(obj.data_type == obj.STRING_ID);
        end
        function read_file(obj, file_name)
            cur_files = dir(file_name);
            if isempty(cur_files)
                error(['Could not find ' ,file_name])
            end
            cur_files([cur_files.isdir]) = [];
            bytes = [cur_files.bytes];
            tot_bytes = sum(bytes);
            obj.raw = zeros(tot_bytes, 1, 'uint8');
            cur_pos = 1;
            for cf = 1 : numel(cur_files)
                obj.files(cf) = string( ...
                    fullfile(cur_files(cf).folder, cur_files(cf).name) ...
                    );
                fid = fopen(obj.files(cf), 'r');
                end_pos = cur_pos + bytes(cf) - 1;
                [obj.raw(cur_pos : end_pos), count] = fread(fid, '*uint8');
                fclose(fid);
                cur_pos = cur_pos + count;
            end
            if cur_pos < tot_bytes
                obj.raw(cur_pos : end) = [];
            end
            obj.find_headers();
        end

        function val = xform(obj,dst,src,varargin)
        end
        function vel = velocity(obj,dst,filter)
        end
    end
    methods(Access=protected)
        function find_headers(obj)
            cur_pos = 0;
            c_ens = 1;
            est_size = floor(obj.nbytes / 1000);
            obj.data_position = nan(est_size, 1);
            obj.data_type = zeros(est_size, 1, 'uint8');
            while cur_pos < obj.nbytes
                cur_pos = cur_pos + 1;
                f_sync = 1;
                if ~obj.raw(cur_pos) == 0xa5
                    f_sync = find(obj.raw(cur_pos : end) == 0xa5, 1, 'first');
                end
                if isempty(f_sync)
                    break
                end
                f_sync = f_sync + cur_pos - 1;
                if f_sync + 1 > obj.nbytes
                    break
                end
                head_size = double(obj.raw(f_sync + 1));
                if f_sync + head_size > obj.nbytes
                    continue
                end
                csum_pos = f_sync + head_size - 2;
                csum_calc = obj.compute_checksum(f_sync, head_size - 2);
                csum_read = typecast( ...
                    obj.raw(csum_pos: csum_pos + 1), 'uint16' ...
                    );
                if csum_calc ~= csum_read
                    continue
                end
                size_data = double(typecast(obj.raw(f_sync + 4 : f_sync + 5), 'uint16'));
                if f_sync + head_size + size_data - 1 > obj.nbytes
                    continue
                end
                csum_calc = obj.compute_checksum(f_sync + head_size, size_data);
                csum_read = typecast(obj.raw(f_sync + 6: f_sync + 7), 'uint16');
                if csum_calc~=csum_read
                    continue
                end
                if c_ens > numel(obj.data_position)
                    obj.data_position(numel(obj.data_position)*2) = nan;
                end
                obj.data_position(c_ens) = f_sync + head_size;
                obj.data_type(c_ens) = obj.raw(f_sync + 2);
                c_ens = c_ens + 1;
                cur_pos = f_sync + head_size + size_data - 1;
            end
            obj.data_position(c_ens:end) = [];
            obj.data_type(c_ens:end) = [];
        end
        function csum = compute_checksum(obj, start, size)
            n_shorts = floor(size/2);
            bytes_left = mod(size,2) == 1;
            csum = double(0xb58c) + sum(...
                double(typecast(...
                obj.raw(start : start + n_shorts * 2 - 1), ...
                'uint16')));
            if bytes_left
                csum = csum + double(obj.raw(start + size - 1));
            end
            csum = uint16(mod(csum, 65536));
        end
        function val = get_nbeams(obj)
        end
        function val = get_nensembles(obj)
            
        end
        function val = get_ncells(obj)
        end
        function val = get_coordinate_system(obj)
        end
        function val = get_beam_angle(obj)
        end
        function val = get_pitch(obj)
        end
        function val = get_roll(obj)
        end
        function val = get_heading(obj)
        end
        function val = get_cellsize(obj)
        end
        function val = get_time(obj)
        end
        function val = get_distmidfirstcell(obj)
        end
        function val = get_temperature(obj)
        end
        function val = get_salinity(obj)
        end
        function val = get_pressure(obj)
        end
        function val = get_echo(obj)
        end
        function val = get_backscatter(obj)
        end
        function val = get_transducer(obj)
        end
    end
    properties (Access = protected)
        raw (1, :) uint8
        data_position(:, 1) double
        data_type(:,1) uint8
    end
    properties (SetAccess = protected)
        files (:, 1) string
    end
end