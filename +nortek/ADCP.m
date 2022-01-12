classdef ADCP < ADCP
    properties(Constant)
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

        BURST_BIT_PRESSURE = 0;
        BURST_BIT_TEMPERATURE = 1;
        BURST_BIT_COMPASS = 2;
        BURST_BIT_TILT = 3;
        BURST_BIT_VELOCITY = 5;
        BURST_BIT_AMPLITUDE = 6;
        BURST_BIT_CORRELATION = 7;
        BURST_BIT_ALTIMETER = 8;
        BURST_BIT_ALTIMETER_RAW = 9;
        BURST_BIT_AST = 10;
        BURST_BIT_ECHO_SOUNDER = 11;
        BURST_BIT_AHRS = 12;
        BURST_BIT_PERC_GOOD = 13;
        BURST_BIT_STDDEV = 14
        
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
        blanking
        velocity_scale
        number_altimeter_samples
        configuration_string
    end
    properties(Dependent, Access = protected)
        burst_version
    end
    methods
        function obj = ADCP(varargin)
            obj = obj@ADCP(varargin{:})
            for ca = 1 : nargin
                if ischar(varargin{ca}) || isStringScalar(varargin{ca})
                    obj.read_file(varargin{ca})
                end
            end
            obj.transformation_matrix_source = nortek.InstrumentMatrixFromConfString;
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
        function val = get.burst_version(obj)
            data_pos = obj.data_position(obj.data_type == obj.BURST_ID);
            val = double(obj.raw(data_pos));
        end
        function val = get.blanking(obj)
            val = double(obj.get_scalar([38 34], 'uint16', obj.BURST_ID) );
            scale = ones(1, obj.nensembles) * 0.01; %TODO: check if this 
                                                    % scaling is correct
            scale(obj.burst_version == 2) = 0.001;
            val = val .* scale;
        end
        function val = get.velocity_scale(obj)
            val = double(obj.get_scalar([62 58], 'int8', obj.BURST_ID) );
        end
        function val = get.number_altimeter_samples(obj)
            if ~burst_has_data(obj.BURST_BIT_ALTIMETER_RAW)
                error('No raw altimeter data available')
            end
            if any(obj.burst_version == 2)
                error('No altimeter data in burst data version 2')
                % maybe here we should get it from outside the burst data
            end
            of = obj.burst_data_offset(obj.BURST_BIT_ALTIMETER_RAW);
            val = double(obj.get_scalar([of of]), 'uint16', obj.BURST_ID);
        end
        function val=get.configuration_string(obj)
            fstring = find(obj.data_type == obj.STRING_ID);
            fstring = fstring(obj.raw(obj.data_position(fstring)) == 0x10);
            if isempty(fstring)
                error('No configuration available')
            end
            fstring=fstring(1);
            pos = obj.data_position(fstring);
            dat_size = obj.data_size(fstring);

            % exclude 0x10 and ending zero byte
            val = char(obj.raw(pos + 1: pos + dat_size - 2)); 
        end
        function val = burst_has_data(obj,bit)
            dat = obj.get_scalar([6 2], 'uint16', obj.BURST_ID);
            val = logical(bitand(1, bitshift(dat, -bit)));
        end
        function offs = burst_data_offset(obj,bit)
            offs = zeros(1, obj.nensembles);
            offs(obj.burst_version == 2) = ...
                offs(obj.burst_version == 2) + 68;
            offs(obj.burst_version == 3) = ...
                offs( obj.burst_version == 3) + 76;
            nc = obj.ncells;
            nfields = obj.nbeams .* nc;
            check_bit = @(x) bit > x & obj.burst_has_data(x);
            offs = offs + check_bit(obj.BURST_BIT_VELOCITY) .* nfields * 2;
            offs = offs + check_bit(obj.BURST_BIT_ECHO_SOUNDER) .* nfields;
            offs = offs + check_bit(obj.BURST_BIT_CORRELATION) .* nfields;           
            offs = offs + check_bit(obj.BURST_BIT_ALTIMETER) * 8;
            offs = offs + check_bit(obj.BURST_BIT_AST) * 20;
            if obj.burst_has_data(obj.BURST_BIT_ALTIMETER_RAW)
                offs = offs + check_bit(obj.BURST_BIT_ALTIMETER_RAW) .*...
                    (4 + 2 * obj.number_altimeter_samples);
            end
            offs = offs + check_bit(obj.BURST_BIT_ECHO_SOUNDER) .* nc * 2;
            offs = offs + check_bit(obj.BURST_BIT_AHRS) * 16 * 4;
            offs = offs + check_bit(obj.BURST_BIT_PERC_GOOD) .* nc;
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
            if ~ obj.burst_has_data(obj.BURST_BIT_VELOCITY)
                error('No velocity data in burst record')
            end
            offs = obj.burst_data_offset(obj.BURST_BIT_VELOCITY);
            vel = double(obj.get_field(offs, 'int16', obj.BURST_ID));
            vel = vel .* 10.^obj.velocity_scale;

            %TODO: implement transformations and filters
        end
    end
    methods(Access=protected)
        function find_headers(obj)
            cur_pos = 0;
            c_ens = 1;
            est_size = floor(obj.nbytes / 1000);
            obj.data_position = nan(est_size, 1);
            obj.data_type = zeros(est_size, 1, 'uint8');
            obj.data_size = nan(est_size, 1);
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
                    obj.data_type(numel(obj.data_type)*2) = nan;
                    obj.data_size(numel(obj.data_size)*2) = nan;
                end
                obj.data_position(c_ens) = f_sync + head_size;
                obj.data_type(c_ens) = obj.raw(f_sync + 2);
                obj.data_size(c_ens) = size_data;
                c_ens = c_ens + 1;
                cur_pos = f_sync + head_size + size_data - 1;
            end
            obj.data_position(c_ens:end) = [];
            obj.data_type(c_ens:end) = [];
            obj.data_size(c_ens:end) = [];
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
        function val = get_scalar(obj, offset, data_type, data_id)
            of = obj.data_position(obj.data_type == data_id);
            of = reshape(of, 1, []);
            ver = obj.burst_version;
            of(ver == 2) = of(ver == 2) + offset(1);
            of(ver == 3) = of(ver == 3) + offset(2);
            of = (0 : obj.sizeof(data_type) -1 )' + of;
            of = of(:);
            val = typecast(obj.raw(of), data_type);
        end
        function val = get_field(obj, offset, data_type, data_id)
            of = obj.data_position(obj.data_type == data_id);
            of = reshape(of, 1, []) + offset;
            siz = obj.sizeof(data_type);
            nc = obj.ncells;
            nb = obj.nbeams;
            nfields = nc .* nb .* siz;
            max_fields = max(nfields);
            nens = obj.nensembles;
            of_add = cumsum(ones(max_fields,nens), 1)-1;
            fgood = of_add < nfields;
            of_add( ~fgood ) = nan;
            of = of + of_add;
            of = of(fgood);
            dat = nan(max(nc), max(nb), nens);
            count_beams = cumsum(ones(max(nc), max(nb), nens), 2)-1;
            count_cells = cumsum(ones(max(nc), max(nb), nens), 1)-1;
            fgood_mat = count_beams < shiftdim(nb, -1) & ...
                        count_cells < shiftdim(nc, -1);
            dat(fgood_mat) = typecast(obj.raw(of), data_type);

            val = reshape(dat, max(nc), max(nb), nens);
            val = permute(val, [1, 3, 2]);
        end
        function val = get_nbeams(obj)
            val = obj.get_scalar([34, 30], 'uint16', obj.BURST_ID);
            val = double(bitshift(val, -12));
        end
        function val = get_nensembles(obj)
            val = sum(obj.data_type==obj.BURST_ID);
        end
        function val = get_ncells(obj)
            val = obj.get_scalar([34, 30], 'uint16', obj.BURST_ID);
            val = double(bitand(val, 2^10-1));
        end
        function val = get_coordinate_system(obj)
            dat = obj.get_scalar([34, 30], 'uint16', obj.BURST_ID);
            dat = double(bitand(bitshift(dat, -10), 2^2-1));
            val(numel(dat))=CoordinateSystem.Beam;
            val(dat == 0) = CoordinateSystem.Earth;
            val(dat == 1) = CoordinateSystem.Instrument;
            val(dat == 2) = CoordinateSystem.Beam;
        end
        function val = get_beam_angle(obj)
            % get it from conf string
        end
        function val = get_heading(obj)
            val = double(obj.get_scalar([24 24], 'uint16', obj.BURST_ID));
            val = val / 100;
            val(~obj.burst_has_data(obj.BURST_BIT_COMPASS)) = nan;
        end
        function val = get_pitch(obj)
            val = double(obj.get_scalar([26 26], 'int16', obj.BURST_ID));
            val = val / 100;
            val(~obj.burst_has_data(obj.BURST_BIT_TILT)) = nan;
        end
        function val = get_roll(obj)
            val = double(obj.get_scalar([28 28], 'int16', obj.BURST_ID));
            val = val / 100;
            val(~obj.burst_has_data(obj.BURST_BIT_TILT)) = nan;
        end
        function val = get_cellsize(obj)
            val = double(obj.get_scalar([36 32], 'uint16', obj.BURST_ID));
            val = val / 1000;
        end
        function val = get_time(obj)
            id = obj.BURST_ID;
            year = double(obj.get_scalar(8 * [1 1], 'uint8', id) );
            year = year + 1900;
            month = double(obj.get_scalar(9 * [1 1], 'uint8', id) ) + 1;
            day = double(obj.get_scalar(10 * [1 1], 'uint8', id) );
            hour = double(obj.get_scalar(11 * [1 1], 'uint8', id) );
            minute = double(obj.get_scalar(12 * [1 1], 'uint8', id) );
            second = double(obj.get_scalar(13 * [1 1], 'uint8', id) );
            msec = double(obj.get_scalar(14 * [1 1], 'uint16', id) );
            msec = msec / 10;
            val = datetime(year, month, day, hour, minute, second, msec);
        end

        function val = get_distmidfirstcell(obj)
            val = obj.cellsize + obj.blanking;
        end
        function val = get_temperature(obj)
            val = double(obj.get_scalar([18 18], 'int16', obj.BURST_ID) );
            val = val / 100;
            val(~obj.burst_has_data(obj.BURST_BIT_TEMPERATURE)) = nan;
        end
        function val = get_salinity(obj)
            val = nan(1, obj.nensembles);
        end
        function val = get_pressure(obj)
            val = double(obj.get_scalar([20 20], 'uint32', obj.BURST_ID) );
            val = val * 10;
            val(~obj.burst_has_data(obj.BURST_BIT_PRESSURE)) = nan;
        end
        function val = get_echo(obj)
            if ~ obj.burst_has_data(obj.BURST_BIT_AMPLITUDE)
                error('Amplitude data not available in burst record')
            end
            offs = obj.burst_data_offset(obj.BURST_BIT_AMPLITUDE);
            val = obj.get_field(offs, 'uint8', obj.BURST_ID);
        end
        function val = get_backscatter(obj)
        end
        function val = get_transducer(obj)
        end
    end
    properties %(Access = protected)
        raw (1, :) uint8
        data_position(:, 1) double
        data_type(:,1) uint8
        data_size(:,1) double
    end
    properties (SetAccess = protected)
        files (:, 1) string
    end
    methods(Static)
        function nb=sizeof(type)
            if strcmp(type,'char')
                nb=1;
                return
            end
            nb=cast(1,type); %#ok<NASGU> 
            nb=whos('nb');
            nb=nb.bytes;
        end
    end
end