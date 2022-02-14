classdef VMADCP < nortek.ADCP & VMADCP
    properties(Dependent, SetAccess = protected)
        gnss_nbytes
        gnss_available_packets
        gnss_nensembles
        gnss_time
        gnss_latitude
        gnss_longitude
        gnss_height
        gnss_vel_east
        gnss_vel_north
        gnss_vel_down
        gnss_acc_x
        gnss_acc_y
        gnss_acc_z
        gnss_roll
        gnss_pitch
        gnss_heading
        gnss_ang_velx
        gnss_ang_vely
        gnss_ang_velz
        gnss_std_latitude
        gnss_std_longitude
        gnss_std_height
        bt_nensembles
        bt_nbeams
        bt_time
        bt_velocity_scale
        bt_velocity
        bt_coordinate_system
        bt_depth
    end
    methods
        function obj = VMADCP(varargin)
            search_string_idx = [];
            search_string = '';
            for ca = 1 : nargin
                if ischar(varargin{ca}) || isStringScalar(varargin{ca})
                    search_string_idx = ca;
                    search_string = varargin{search_string_idx};
                end
            end
            varargin{search_string_idx} = [];
            obj = obj@nortek.ADCP(varargin{:})
            obj = obj@VMADCP(varargin{:})
            if ~isempty(search_string)
                obj.read_sigvm(search_string)
            end
        end
        function val = get.gnss_nensembles(obj)
            val = sum(obj.gnss_packet_id == nortek.AnppPacketId.SystemState);
        end
        function val = get.gnss_time(obj)
            unix_time = double(get_gnss_scalar(obj, 4, 'uint32'))+...
                double(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 'uint32'))/1e6;
            val = datetime(unix_time,'ConvertFrom','posixtime','TimeZone','UTC'); % this goes wrong with leap seconds!
        end
        function val = get.gnss_latitude(obj)
            val = rad2deg(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 12, 'double'));
        end
        function val = get.gnss_longitude(obj)
            val = rad2deg(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 20, 'double'));
        end
        function val = get.gnss_height(obj)
            val = get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 28, 'double');
        end
        function val = get.gnss_vel_east(obj)
            val = double(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 40, 'single'));
        end
        function val = get.gnss_vel_north(obj)
            val = double(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 36, 'single'));
        end
        function val = get.gnss_vel_down(obj)
            val = double(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 44, 'single'));
        end
        function val = get.gnss_acc_x(obj)
            val = double(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 48, 'single'));
        end
        function val = get.gnss_acc_y(obj)
            val = double(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 52, 'single'));
        end
        function val = get.gnss_acc_z(obj)
            val = double(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 56, 'single'));
        end
        function val = get.gnss_roll(obj)
            val = rad2deg(double(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 64, 'single')));
        end
        function val = get.gnss_pitch(obj)
            val = rad2deg(double(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 68, 'single')));
        end
        function val = get.gnss_heading(obj)
            val = rad2deg(double(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 72, 'single')));
        end
        function val = get.gnss_ang_velx(obj)
            val = rad2deg(double(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 76, 'single')));
        end
        function val = get.gnss_ang_vely(obj)
            val = rad2deg(double(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 80, 'single')));
        end
        function val = get.gnss_ang_velz(obj)
            val = rad2deg(double(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 84, 'single')));
        end
        function val = get.gnss_std_latitude(obj)
            val = double(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 88, 'single'));
        end
        function val = get.gnss_std_longitude(obj)
            val = double(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 92, 'single'));
        end
        function val = get.gnss_std_height(obj)
            val = double(get_gnss_scalar(obj,...
                nortek.AnppPacketId.SystemState, 96,'single'));
        end
        function val = get.gnss_nbytes(obj)
            val = numel(obj.gnss_raw);
        end
        function val = get.gnss_available_packets(obj)
            val = unique(obj.gnss_packet_id);
        end
        function val = get_bt_filt(obj)
            val = get_data_filt(obj, nortek.DataHeader.BottomTracking,...
                obj.configuration_selector);
        end
        function val = get.bt_nensembles(obj)
            filt = obj.get_bt_filt;
            val = obj.get_nensembles(filt);
        end
        function val = get.bt_nbeams(obj)
            filt = obj.get_bt_filt;
            val = get_nbeams(obj,filt);
        end

        function val = get.bt_time(obj)
            filt = obj.get_bt_filt;
            val = obj.get_time(filt);
        end

        function val = get.bt_velocity_scale(obj)
            filt = obj.get_bt_filt;
            val = double(obj.get_scalar(60, 'int8', filt));
        end

        function val = get.bt_velocity(obj)
             val = read_btfield(obj, 78, 'int32');
             fbad = val == int32(10.^(-obj.bt_velocity_scale)*-9.9);
             val = double(val) .* 10 .^(obj.bt_velocity_scale);
             val(fbad) = nan;
        end

        function val = get.bt_coordinate_system(obj)

        end

        function val = get.bt_depth(obj)

        end

    end

    methods(Access = protected)
        function val = read_btfield(obj,offset, data_type)
            filt = obj.get_bt_filt;
            of = obj.data_position(filt);
            of = reshape(of, 1, []) + offset;
            siz = obj.sizeof(data_type);
            nb = obj.get_nbeams(filt);
            nc = ones(size(nb));
            nfields = nc .* nb .* siz;
            max_fields = max(nfields);
            nens = obj.get_nensembles(filt);
            of_add = cumsum(ones(max_fields,nens), 1)-1;
            of = of + of_add;
            of = of(of_add < nfields);
            clearvars of_add
            count_beams = cumsum(ones(max(nc), max(nb), nens), 2)-1;
            count_cells = cumsum(ones(max(nc), max(nb), nens), 1)-1;
            fgood_mat = count_beams < shiftdim(nb, -1) & ...
                        count_cells < shiftdim(nc, -1);
            clearvars count_beams count_cells
            val = nan(max(nc), max(nb), nens);
            val(fgood_mat) = typecast(obj.raw(of), data_type);
            val = permute(val, [1, 3, 2]);
        end

        function read_sigvm(obj, search_string)
            cur_files = dir(search_string);
            if isempty(cur_files)
                error(['Could not find ' ,search_string])
            end
            cur_files([cur_files.isdir]) = [];

            unzip_dest = tempname;
            if exist(unzip_dest,"dir")
                error('cannot create temporary directory');
            else
                mkdir(unzip_dest)
            end
            nfiles = numel(cur_files);

            ad2cp_files = cell(nfiles,1);
            gnss_files = cell(nfiles,1);
            for cf = 1:numel(cur_files)
                try
                    dst=unzip(fullfile(cur_files(cf).folder, cur_files(cf).name), unzip_dest);
                catch err
                    error([cur_files(cf).name,' is not a valid SigVM file.'])
                end
                % store gps and ad2cp file names and read them.
                fgps = find(contains(dst, '.anpp','IgnoreCase',true), 1);
                gnss_files(cf) = dst(fgps);
                if isempty(fgps)
                    error(['No gps data in file:', cur_files(cf).name])
                end
                fad2cp = find(contains(dst, '.ad2cp','IgnoreCase',true), 1);
                if isempty(fad2cp)
                    error(['No adcp data in file:', cur_files(cf).name])
                end
                ad2cp_files(cf) = dst(fad2cp);
            end
            disp('reading gnss files...')
            obj.read_gnss(gnss_files);
            disp('reading ad2cp files...')
            obj.read_ad2cp(ad2cp_files);


            % cleanup temporary files
            disp('cleaning up temporary files...')
            delete(fullfile(unzip_dest, '*'))
            rmdir(unzip_dest)

            disp('done.')
        end
        function read_gnss(obj, files)
            obj.gnss_raw = obj.read_file(files);
            pos = 1;
            est_size = floor(obj.gnss_nbytes / 32);
            store_pos = zeros(est_size, 1, 'uint32');
            store_id = zeros(est_size, 1, 'uint8');
            store_size = zeros(est_size, 1,'uint8');
            c_ens=1;
            tot_bytes = obj.gnss_nbytes;
            while pos + 4 <= numel(obj.gnss_raw)
                if ~obj.valid_lrc(pos)
                    disp(['bad lrc at position ', num2str(pos)])
                    pos = pos + 1;
                    continue
                end
                % check buffer size is large enough for crc computation
                tmp_size = obj.gnss_raw(pos+2);
                if pos + 4 + double(tmp_size) > tot_bytes ||...
                        ~obj.valid_crc(pos, tmp_size)
                    disp(['bad crc at position ', num2str(pos)])
                    pos = pos + 1;
                    continue
                end
                % valid packet is found!

                % expand storage variables if necessary
                if c_ens > est_size
                    disp('expanding')
                    est_size = est_size * 2;
                    store_pos(est_size) = 0x00000000;
                    store_id(est_size) = 0xff;
                    store_size(est_size) = 0x00;
                end

                % store packet information
                store_pos(c_ens) = pos+5;
                store_id(c_ens) = obj.gnss_raw(pos+1);
                store_size(c_ens) = tmp_size;
                c_ens = c_ens + 1;

                % move position to start of new data packet
                pos = pos + 5 + double(tmp_size);
            end

            % trim output variables
            store_pos(c_ens:end) = [];
            store_id(c_ens:end) = [];
            store_size(c_ens:end) = [];

            % store outputs in object
            obj.gnss_position = store_pos;
            obj.gnss_packet_id = store_id;
            obj.gnss_size = store_size;
        end
        function val = valid_lrc(obj,pos)
            val = obj.gnss_raw(pos) == uint8(bitand(bitxor(bitand(sum(uint16(obj.gnss_raw(pos+(1:4))),'native'),0x00ff),0x00ff)+1,0x00ff));
        end
        function val = valid_crc(obj, pos, packet_size)
            crc = 0xffff;
            for cb = pos + 5 + (0: double(packet_size) - 1)
                idx = bitxor(obj.gnss_raw(cb),...
                    uint8(bitshift(crc, -8)));
                crc = bitxor(obj.crc_table(double(idx) + 1),...
                    bitshift(crc, 8));
            end
            val = crc == typecast(obj.gnss_raw(pos + (3:4)),'uint16');
        end
        function val = get_gnss_scalar(obj, packet, offs, type)
            packet_filt = obj.gnss_packet_id == packet;
            of = obj.gnss_position(packet_filt)';
            of = offs + uint32((0 : obj.sizeof(type) -1 ))' + of;
            of = of(:);
            val = typecast(obj.gnss_raw(of), type);
        end
        function val = get_btvel(obj)
        end
        function val = get_bt_vertical_range(obj)	% defined in VMADCP
        end
    end
    properties%(Access = protected)
        gnss_raw (1, :) uint8
        gnss_packet_id (:,1) nortek.AnppPacketId
        gnss_size (:,1) uint8
        gnss_position (:,1) uint32
    end
    properties(Constant, Access = protected)
        crc_table = [
            0x0000, 0x1021, 0x2042, 0x3063, 0x4084, 0x50a5,...
            0x60c6, 0x70e7, 0x8108, 0x9129, 0xa14a, 0xb16b,...
            0xc18c, 0xd1ad, 0xe1ce, 0xf1ef, 0x1231, 0x0210,...
            0x3273, 0x2252, 0x52b5, 0x4294, 0x72f7, 0x62d6,...
            0x9339, 0x8318, 0xb37b, 0xa35a, 0xd3bd, 0xc39c,...
            0xf3ff, 0xe3de, 0x2462, 0x3443, 0x0420, 0x1401,...
            0x64e6, 0x74c7, 0x44a4, 0x5485, 0xa56a, 0xb54b,...
            0x8528, 0x9509, 0xe5ee, 0xf5cf, 0xc5ac, 0xd58d,...
            0x3653, 0x2672, 0x1611, 0x0630, 0x76d7, 0x66f6,...
            0x5695, 0x46b4, 0xb75b, 0xa77a, 0x9719, 0x8738,...
            0xf7df, 0xe7fe, 0xd79d, 0xc7bc, 0x48c4, 0x58e5,...
            0x6886, 0x78a7, 0x0840, 0x1861, 0x2802, 0x3823,...
            0xc9cc, 0xd9ed, 0xe98e, 0xf9af, 0x8948, 0x9969,...
            0xa90a, 0xb92b, 0x5af5, 0x4ad4, 0x7ab7, 0x6a96,...
            0x1a71, 0x0a50, 0x3a33, 0x2a12, 0xdbfd, 0xcbdc,...
            0xfbbf, 0xeb9e, 0x9b79, 0x8b58, 0xbb3b, 0xab1a,...
            0x6ca6, 0x7c87, 0x4ce4, 0x5cc5, 0x2c22, 0x3c03,...
            0x0c60, 0x1c41, 0xedae, 0xfd8f, 0xcdec, 0xddcd,...
            0xad2a, 0xbd0b, 0x8d68, 0x9d49, 0x7e97, 0x6eb6,...
            0x5ed5, 0x4ef4, 0x3e13, 0x2e32, 0x1e51, 0x0e70,...
            0xff9f, 0xefbe, 0xdfdd, 0xcffc, 0xbf1b, 0xaf3a,...
            0x9f59, 0x8f78, 0x9188, 0x81a9, 0xb1ca, 0xa1eb,...
            0xd10c, 0xc12d, 0xf14e, 0xe16f, 0x1080, 0x00a1,...
            0x30c2, 0x20e3, 0x5004, 0x4025, 0x7046, 0x6067,...
            0x83b9, 0x9398, 0xa3fb, 0xb3da, 0xc33d, 0xd31c,...
            0xe37f, 0xf35e, 0x02b1, 0x1290, 0x22f3, 0x32d2,...
            0x4235, 0x5214, 0x6277, 0x7256, 0xb5ea, 0xa5cb,...
            0x95a8, 0x8589, 0xf56e, 0xe54f, 0xd52c, 0xc50d,...
            0x34e2, 0x24c3, 0x14a0, 0x0481, 0x7466, 0x6447,...
            0x5424, 0x4405, 0xa7db, 0xb7fa, 0x8799, 0x97b8,...
            0xe75f, 0xf77e, 0xc71d, 0xd73c, 0x26d3, 0x36f2,...
            0x0691, 0x16b0, 0x6657, 0x7676, 0x4615, 0x5634,...
            0xd94c, 0xc96d, 0xf90e, 0xe92f, 0x99c8, 0x89e9,...
            0xb98a, 0xa9ab, 0x5844, 0x4865, 0x7806, 0x6827,...
            0x18c0, 0x08e1, 0x3882, 0x28a3, 0xcb7d, 0xdb5c,...
            0xeb3f, 0xfb1e, 0x8bf9, 0x9bd8, 0xabbb, 0xbb9a,...
            0x4a75, 0x5a54, 0x6a37, 0x7a16, 0x0af1, 0x1ad0,...
            0x2ab3, 0x3a92, 0xfd2e, 0xed0f, 0xdd6c, 0xcd4d,...
            0xbdaa, 0xad8b, 0x9de8, 0x8dc9, 0x7c26, 0x6c07,...
            0x5c64, 0x4c45, 0x3ca2, 0x2c83, 0x1ce0, 0x0cc1,...
            0xef1f, 0xff3e, 0xcf5d, 0xdf7c, 0xaf9b, 0xbfba,...
            0x8fd9, 0x9ff8, 0x6e17, 0x7e36, 0x4e55, 0x5e74,...
            0x2e93, 0x3eb2, 0x0ed1, 0x1ef0];
    end
end