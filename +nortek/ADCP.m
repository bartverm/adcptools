classdef ADCP < ADCP
% Nortek ADCP dataset
%   
%   nortek.ADCP is an ADCP class that represents Nortek ADCP data.
%
%   A = nortek.ADCP('search string') opens ad2cp files found with the
%   search string passed to the object
%
%   nortek.ADCP properties (configuration):
%   data_selector - select data type(s) to be used
%   configuration_selector - select configuration(s) to be used
%   echo_sounder_selector - select sounder configuration(s) to be used
%
%   ADCP properties (configuration, inherited):
%   filters - filters to exclude data from profiles
%   timezone - specifies the timezone used
%   horizontal_position_provider - provides horizontal positioning
%   vertical_position_provider - provides vertical positioning
%   tilts_provider - provides pitch and roll angles
%   heading_provider - provides heading angle
%   instrument_matrix_provider - provides instrument to beam matrix
%   transducer - transducer acoustic properties
%   water - water acoustic properties
%
%   nortek.ADCP properties (read only):
%   available_data - types of data available in the file
%   available_configurations - configurations available in the file
%   available_echo_sounder - echo sounder configurations available
%   nbytes - number of bytes in raw data
%   has_heading_internal - whether internal heading is available
%   heading_internal - internally recorder instrument heading
%   has_tilts_internal - whether internal pitch and roll are available
%   pitch_internal - internally recorder pitch angle
%   roll_internal - internally recorder roll angle
%   blanking - blanking distance (m)
%   number_altimeter_samples - number of altimeter samples
%   configuration_string - configuration during data collection
%   ahrs_matrix - ahrs orientation matrix
%   sounder_ncells - number of cells used by the echo sounder
%   sounder_nensembles - number of ensembles used by echo sounder
%   sounder_time - time of echo sounder ensembles
%   sounder_cellsize - cell size of echo sounder cells (m)
%   sounder_nbeams - number of beams used by echo sounder
%   sounder_echo - received echo strength of echo sounder (dB)
%
%   ADCP data properties (read only, inherited):
%   nbeams - number of acoustic beams
%   nensembles - number of ensembles
%   ncells - number of depth cells
%   coordinate_system - coordinate system of velocity data
%   beam_angle - angle of acoustic beams with vertical (degrees)
%   pitch - pitch angle (degrees)
%   roll - roll angle (degrees)
%   heading - heading angle (degrees)
%   cellsize - depth cell size (m)
%   time - time and date of ensembles
%   distmidfirstcell - distance to center of first depth cell (m)
%   depth_cell_slant_range - distance along ac. beam to depth cells (m)
%   temperature - instrument temperature (Celsius)
%   salinity - salinity of water (psu)
%   pressure - pressure of water (Pa)
%   echo - received echo intensity (dB)
%   backscatter - volume backscatter strength (dB)
%   horizontal_position - horizontal position of the instrument (m)
%   vertical_position - vertical position of the instrument (m)
%   beam_2_instrument_matrix - beam to instrument transformation matrix
%   instrument_2_beam_matrix - instrument to beam transformation matrix
%
%   nortek.ADCP methods:
%   read_file - read content of an ad2cp file
%
%   ADCP methods (inherited):
%   bad - return a mask with ones for data marked as being bad
%   plot_orientations - plots pitch, roll and heading of the ADCP
%   plot_velocity - plot the velocity in Earth coordinates
%   plot_backscatter - plot the volume backscatter strength
%   plot_all - plot all available plots
%   depth_cell_offset - get east,north,up offset of depth cells from ADCP
%   depth_cell_position - get east,north,up position of depth cells
%   xform - get transformation matrices
%   velocity - get velocity
%
%   see also: ADCP, rdi.ADCP

    properties
        % Select data to be included. This can be a column vector with
        % different values to include different types of data. Values must 
        % be convertible to nortek.DataHeader type.
        %
        %  see also: nortek.ADCP, available_data, nortek.DataHeader
        data_selector (:,1) nortek.DataHeader = nortek.DataHeader.Burst

        % Select configuration to be used. This can be a column vector with
        % different configurations to be used. 
        %
        %  see also: nortek.ADCP, available_configurations
        configuration_selector (:,1) double = 0

        % Select the echo sounder configuration to use. This can be a
        % column vector with different values to combine different echo
        % sounder configurations
        %
        %  see also: nortek.ADCP, available_echo_sounder
        echo_sounder_selector (:,1) double = 0

        head_tilt_provider (:,1) nortek.HeadTiltProvider = nortek.HeadTiltFromAHRS
    end
    properties(Dependent, SetAccess = private)
        % Available data types in ad2cp file. Column vector of type:
        % nortek.DataHeader
        %
        %   see also: nortek.ADCP, data_selector
        available_data (:,1) nortek.DataHeader

        % Available configurations in the ad2cp file. This is usefull when
        % different configurations are used in the ad2cp file.
        %
        %   see also: nortek.ADCP, configuration_selector
        available_configurations (:,1) double

        % Available echo sounder settings. Echo sounder can be used with
        % different settings. This returns the available conifguragions
        % used.
        %
        %   see also: nortek.ADCP, echo_sounder_selector
        available_echo_sounder (:,1) double

        % Number of bytes in raw data.
        %
        %   see also: nortek.ADCP
        nbytes

        % Whether internally measured heading is available
        %
        %   see also: nortek.ADCP, heading_internal, heading
        has_heading_internal

        % Azimuth angle in degrees of the instrument heading
        %
        % see also: nortek.ADCP, has_heading_internal
        heading_internal

        % Whether internally measured pitch and roll are available.
        %
        % see also: nortek.ADCP, roll_intenal, pitch_internal
        has_tilts_internal

        % Internally measured pitch angle in degrees.
        %
        %   see also: nortek.ADCP, has_tilts_internal
        pitch_internal

        % Internally measured roll angle in degrees.
        %
        %   see also: nortek.ADCP, has_tilts_internal
        roll_internal

        
        head_tilt


        % Blanking distance (m).
        %
        %   see also: nortek.ADCP, blanking_scaling
        blanking

        % Number of altimeter samples.
        %
        %   see also: nortek.ADCP
        number_altimeter_samples

        % Configuration string holding all settings of the instrument
        % during data collection
        %
        %   see also: nortek.ADCP
        configuration_string

        physical_beams_used

        beam_azimuth

        % Rotation matrix based on the AHRS sensor.
        %
        %   see also: nortek.ADCP
        ahrs_matrix

        % Number of cells used by the echo sounder.
        %
        %   see also: nortek.ADCP, sounder_echo
        sounder_ncells

        % Number of echo sounder ensembles
        %
        %   see also: nortek.ADCP, sounder_echo
        sounder_nensembles

        % Time of echo_sounder ensembles
        %
        %   see also: nortek.ADCP, sounder_echo
        sounder_time

        % Cell size of echo sounder (m)
        %
        %   see also: nortek.ADCP, sounder_echo
        sounder_cellsize

        % Number of beams used by sounder
        %
        %   see also: nortek.ADCP, sounder_echo
        sounder_nbeams

        % Received echo strength of echo sounder (dB)
        %
        %   see also: nortek.ADCP
        sounder_echo

        altimeter_distance
    end
    properties(Dependent, SetAccess = private, GetAccess = protected)
        active_configuration (:,1) double
        echo_sounder_index (:,1) double
        blanking_scale
        velocity_scale
        burst_version
    end
    methods
        %%% CONSTRUCTOR
        function obj = ADCP(varargin)
            obj = obj@ADCP(varargin{:})
            for ca = 1 : nargin
                if ischar(varargin{ca}) || isStringScalar(varargin{ca})
                    obj.read_ad2cp(varargin{ca})
                end
            end
            obj.instrument_matrix_provider = [...
                nortek.InstrumentMatrixFromConfString;...
                nortek.InstrumentMatrixFromBAngle];
            obj.heading_provider = nortek.HeadingInternal;
            obj.tilts_provider = nortek.TiltsInternal;
            obj.timezone='UTC';
        end

        %%% SET AND GET METHODS
        function val = get.available_data(obj)
            val = unique(obj.data_header);
        end
        function val = get.available_configurations(obj)
            val = unique(obj.active_configuration);
            val(isnan(val)) = [];
        end
        function val = get.available_echo_sounder(obj)
            val = unique(obj.echo_sounder_index);
            val(isnan(val)) = [];
        end
        function val = get.nbytes(obj)
            val = numel(obj.raw);
        end
        function val = get.burst_version(obj)
            val = obj.get_burst_version;
        end
        function val = get.has_heading_internal(obj)
            val = obj.get_has_heading_internal;
        end
        function val = get.heading_internal(obj)
            val = obj.get_heading_internal;
        end
        function val=get.has_tilts_internal(obj)
            val = obj.get_has_tilts_internal;
        end
        function val = get.pitch_internal(obj)
            val = obj.get_pitch_internal;
        end
        function val = get.roll_internal(obj)
            val = obj.get_roll_internal;
        end
        function val = get.head_tilt(obj)
            val = obj.head_tilt_provider.head_tilt_matrix(obj);
        end
        function val = get.blanking(obj)
            val = obj.get_blanking;
        end
        function val = get.blanking_scale(obj)
            val = obj.get_blanking_scaling;
        end

        function val = get.velocity_scale(obj)
            val = obj.get_velocity_scale;
        end

        function val = get.number_altimeter_samples(obj)
            val = obj.get_number_altimeter_samples;
        end

        function val=get.configuration_string(obj)
            fstring = find(obj.data_header == nortek.DataHeader.String);
            fstring = fstring(obj.raw(obj.data_position(fstring)) == ...
                nortek.DataHeader.ConfigurationString);
            if isempty(fstring)
                error('No configuration available')
            end
            fstring=fstring(1);
            pos = obj.data_position(fstring);
            dat_size = obj.data_size(fstring);

            % exclude 0x10 and ending zero byte
            val = char(obj.raw(pos + 1: pos + dat_size - 2)); 
        end
        function val = get.physical_beams_used(obj)
            val = obj.get_physical_beams_used;
        end
        function val = get.beam_azimuth(obj)
            val = obj.get_beam_azimuth;
        end
        function val = get.ahrs_matrix(obj)
            val = obj.get_ahrs_matrix;
        end

        function val = get.sounder_time(obj)
            val = obj.get_sounder_time;
        end

        function val = get.sounder_nensembles(obj)
            val = obj.get_sounder_nensembles;
        end
        function val = get_sounder_nensembles(obj, filt)
            if nargin < 2
                filt = obj.get_sounder_filt;
            end
            val = obj.get_nensembles(filt);
        end
        function val = get.sounder_cellsize(obj)
            val = obj.get_sounder_cellsize;
        end

        function val = get.sounder_ncells(obj)
            val = obj.get_sounder_ncells;
        end

        
        function val = get.sounder_echo(obj)
            val = obj.get_sounder_echo;
        end

        function val = get.sounder_nbeams(obj)
            val = ones(1,obj.sounder_nensembles);
        end

        function val = get.altimeter_distance(obj)
            val = obj.get_altimeter_distance;
        end
       
        function val = get.active_configuration(obj)
            % only computed for version 3 burst or average data
            val = nan(size(obj.data_header));
            filt = obj.data_header == nortek.DataHeader.Burst |...
                   obj.data_header == nortek.DataHeader.Average |...
                   obj.data_header == nortek.DataHeader.BottomTracking;
            filt_ver = false(size(filt));
            ver = obj.get_burst_version(filt);
            filt_ver(filt) = ver == 3 | ver == 1; % 1 for bottom tracking
            tmp_val = get_scalar(obj, 68, 'uint32', filt_ver);
            val(filt_ver) = double(obj.get_bit(tmp_val, 16));
        end
        function val = get.echo_sounder_index(obj)
            val = nan(size(obj.data_header));
            filt = obj.data_header ~= nortek.DataHeader.String;
            tmp_val = get_scalar(obj,68, 'uint32', filt);
            val(filt) = double(obj.get_bit(tmp_val, 12, 15));
        end

        %%% ORDINARY METHODS
        function read_ad2cp(obj, file_name)
            obj.raw = obj.read_file(file_name);
            obj.find_headers();
        end
        
        %%% INHERITED ABSTRACT METHODS
        function tm = xform(obj,dst,src,varargin)
            if nargin < 3
                src=obj.coordinate_system;
            end
            tm = repmat(shiftdim(eye(4),-2),1,obj.nensembles);
            I=CoordinateSystem.Instrument;
            E=CoordinateSystem.Earth;
            B=CoordinateSystem.Beam;
            exp_cfilt=true(1,obj.nensembles); % dummy filter to expand scalar input
            
            assert(all(dst~=CoordinateSystem.Ship),'Ship coordinate system not supported for nortek ADCPs')
            assert(all(src~=CoordinateSystem.Ship),'Ship coordinate system not supported for nortek ADCPs')

%             assert(all(src >= dst) | all(src<=dst),'It is not possible to combine forward and backward transformations')

            % FORWARD
            % from lower than instrument to instrument
            cfilt = exp_cfilt & dst >= I & src < I;
            tmptm=obj.beam_2_instrument_matrix(:,cfilt,:,:);
            tm(:,cfilt,:,:)=tmptm(:,cfilt,:,:);
            
            % from lower than earth to earth
            cfilt = exp_cfilt & dst >= E & src < E;
            tm(:,cfilt,:,:)=helpers.matmult(...
                obj.head_tilt(:,cfilt,:,:),...
                tm(:,cfilt,:,:));
            
            % INVERSE           
            % from higher than instrument to instrument
            cfilt = exp_cfilt & dst <= I & src > I;               
            tm(1,cfilt,:,:)=helpers.matmult(...
                permute(obj.head_tilt(:,cfilt,:,:),[1,2,4,3]),...
                tm(1,cfilt,:,:));
            
            % from higher than beam to beam
            cfilt = exp_cfilt & dst == B & src > B;
            tm(:,cfilt,:,:)=helpers.matmult(...
                obj.instrument_2_beam_matrix(:,cfilt,:,:),...
                tm(1,cfilt,:,:));
        end
        function vel = velocity(obj,dst,filter,filt)
            if nargin < 4
                filt = obj.get_data_filt;
            end
            if ~ obj.burst_has_data(nortek.BurstBit.Velocity, filt)
                error('No velocity data in burst record')
            end
            offs = obj.burst_data_offset(nortek.BurstBit.Velocity, filt);
            vel = obj.get_field(offs, 'int16', filt);
            sc = obj.get_velocity_scale(filt);
            fbad = vel == int32(10.^(-sc)*-9.9); % -9.9 is bad
            vel = double(vel) .* 10 .^sc;
            vel(fbad) = nan;

            % coordinate transformation
            src = obj.coordinate_system;
            if nargin > 1 && ~all(dst == src)
                tm = obj.xform(dst);
                vel = helpers.matmult(tm, vel);
            end
            
            % filtering
            if nargin > 2
                bad=obj.bad(filter);
            else
                bad=obj.bad();
            end
            vel(bad)=nan;
        end
        function val = burst_has_data(obj, bit, filt)
            if nargin < 3
                filt = obj.get_data_filt;
            end
            dat = obj.get_scalar([6 2], 'uint16', filt);
            val = logical(obj.get_bit(dat, bit));
        end
    end
    methods(Access=protected)
        %%% OVERLOADABLE SET/GET METHODS
        function val = get_blanking(obj,filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            val = double(obj.get_scalar([38 34], 'uint16', filt));
            val = val .* obj.get_blanking_scaling(filt);
        end
        function val = get_blanking_scaling(obj, filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            val = 1e-3 * ones(1,sum(filt));
            stat_bytes = obj.get_scalar([64 68], 'uint32', filt);
            is_cm = logical(obj.get_bit(stat_bytes,1));
            val(obj.get_burst_version == 3 & is_cm) = 1e-2; 
        end
        function val = get_velocity_scale(obj, filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            val = double(obj.get_scalar([62 58], 'int8', filt));
        end
        function val = get_number_altimeter_samples(obj,filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            val = zeros(1, sum(filt));
            if ~obj.burst_has_data(nortek.BurstBit.AltimeterRaw, filt)
                return
            end
            if any(obj.get_burst_version(filt) == 2)
                return
                % maybe here we should get it from outside the burst data
            end
            of = obj.burst_data_offset(nortek.BurstBit.AltimeterRaw, filt);
            val = double(obj.get_scalar([of of], 'uint16', filt));
        end
        function val = get_physical_beams_used(obj, filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            val = nan(1,obj.get_nensembles(filt), 5);
            dat = get_scalar(obj, [58 54], 'uint16',filt);
            ver = obj.get_burst_version;
            v2 = ver == 2;
            v3 = ver == 3;
            for cb = 1:4
                val(1, v2, cb) = double(obj.get_bit(...
                    dat(v2),...
                    (cb - 1) * 3,...
                    cb * 3 -1 ) );
                val(1, v3, cb) =double(obj.get_bit(...
                    dat(v3),...
                    (cb - 1) * 4,...
                    cb * 4 -1 ) );
            end
            cb = 5;
                val(1, v2, cb) = double(obj.get_bit(...
                    dat(v2),...
                    (cb - 1) * 3,...
                    cb * 3 -1 ) );
            val(:, :, all(isnan(val) | val == 0, 2)) = [];
        end
        function val = get_beam_azimuth(obj,filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            conf = obj.configuration_string;
            lines=strsplit(conf,{'\r','\n'});
            f_lines = ~cellfun(@isempty,regexp(lines,'BEAMCFGLIST'));
            lines = lines(f_lines);
            val = regexp(lines,'PHI=(-?\d+\.\d+)', 'tokens');
            val = [val{:}];
            val = cellfun(@str2double,val);
            beam = regexp(lines,'HWBEAM=(\d)', 'tokens');
            beam = [beam{:}];
            beam = cellfun(@str2double,beam);
            [~, beam_idx] = sort(beam);
            val = val(beam_idx(obj.get_physical_beams_used(filt)));
        end
        function val = get_ahrs_matrix(obj, filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            if ~ obj.burst_has_data(nortek.BurstBit.AHRS, filt)
                error('No AHRS data are available')
            end
            offs = obj.burst_data_offset(nortek.BurstBit.AHRS, filt);
            offs = offs + reshape(obj.data_position(filt),1,[]);
            offs = offs + (0:9*4-1)';
            offs = offs(:);
            val = typecast(obj.raw(offs), 'single');
            val = reshape(val, 3, 3, obj.nensembles);
            val = permute(val, [4, 3, 2, 1]);
            val = double(val);
        end
        function val=get_sounder_time(obj, filt)
            if nargin < 2
                filt = obj.get_sounder_filt;
            end
            val = obj.get_time(filt);
        end
        function val = get_sounder_cellsize(obj,filt)
            if nargin < 2
                filt = obj.get_sounder_filt;
            end
            val = obj.get_cellsize(filt);
        end
        function val = get_sounder_ncells(obj,filt)
            if nargin < 2
                filt = obj.get_sounder_filt;
            end
            val = obj.get_ncells(filt);
            val(obj.data_header(filt) ~= ...
                nortek.DataHeader.EchoSounder) = 0;
        end     
        function val = get_sounder_echo(obj,filt)
            if nargin < 2
                filt = obj.get_sounder_filt;
            end
            if ~ obj.burst_has_data(nortek.BurstBit.EchoSounder, filt)
                error('Echo sounder data not available')
            end
            offs = obj.burst_data_offset(nortek.BurstBit.EchoSounder, filt);
            val = double(obj.get_field(offs, 'uint16', filt))*0.01;
        end
        function val = get_heading_internal(obj, filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            val = double(obj.get_scalar([24 24], 'uint16', filt));
            val = val / 100;
            val(~obj.burst_has_data(nortek.BurstBit.Compass, filt)) = nan;
        end
        function val = get_has_heading_internal(obj,filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            val = any(obj.burst_has_data(nortek.BurstBit.Compass,filt));
        end
        function val = get_pitch_internal(obj, filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            val = double(obj.get_scalar([26 26], 'int16', filt));
            val = val / 100;
            val(~obj.burst_has_data(nortek.BurstBit.Tilt, filt)) = nan;
        end
        function val = get_has_tilts_internal(obj,filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            val = any(obj.burst_has_data(nortek.BurstBit.Tilt,filt));
        end
        function val = get_roll_internal(obj, filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            val = double(obj.get_scalar([28 28], 'int16', filt));
            val = val / 100;
            val(~obj.burst_has_data(nortek.BurstBit.Tilt, filt)) = nan;
        end
        function val = get_burst_version(obj,filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            data_pos = obj.data_position(filt);
            val = double(obj.raw(data_pos));
        end
        function val = get_altimeter_distance(obj,filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            of = obj.burst_data_offset(nortek.BurstBit.Altimeter, filt);
            val = double(get_scalar(obj,of,'single',filt));
        end
        

        %%% INHERITED ABSTRACT GET/SET METHODS
        function val = get_nbeams(obj, filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            val = obj.get_scalar([34, 30], 'uint16', filt);
            val = double(obj.get_bit(val, 12, 15));
            val(obj.data_header(filt) == nortek.DataHeader.EchoSounder) = 1;
        end
        function val = get_nensembles(obj, filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            val = sum(filt);
        end
        function val = get_ncells(obj, filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            tmp_val = obj.get_scalar([34, 30], 'uint16', filt);
            val = nan(1, sum(filt));
            sounder_filt = obj.data_header(filt) == nortek.DataHeader.EchoSounder;
            val(sounder_filt) = double(tmp_val(sounder_filt));
            val(~sounder_filt) = obj.get_bit(tmp_val(~sounder_filt),0, 9);
        end
        function val = get_coordinate_system(obj, filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            dat = obj.get_scalar([34, 30], 'uint16', filt);
            dat = double(obj.get_bit(dat, 10, 11));
            val(numel(dat))=CoordinateSystem.Beam;
            val(dat == 0) = CoordinateSystem.Earth;
            val(dat == 1) = CoordinateSystem.Instrument;
            val(dat == 2) = CoordinateSystem.Beam;
        end
        function val = get_beam_angle(obj,filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            conf = obj.configuration_string;
            lines=strsplit(conf,{'\r','\n'});
            f_lines = ~cellfun(@isempty,regexp(lines,'BEAMCFGLIST'));
            lines = lines(f_lines);
            val = regexp(lines,'THETA=(-?\d+\.\d+)', 'tokens');
            val = [val{:}];
            val = cellfun(@str2double,val);
            beam = regexp(lines,'HWBEAM=(\d)', 'tokens');
            beam = [beam{:}];
            beam = cellfun(@str2double,beam);
            [~, beam_idx] = sort(beam);
            val = val(beam_idx(obj.get_physical_beams_used(filt)));
        end

        function val = get_cellsize(obj,filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            val = double(obj.get_scalar([36 32], 'uint16', filt));
            val = val / 1000;
        end
        function val = get_time(obj,filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            year = double(obj.get_scalar(8 * [1 1], 'uint8', filt) );
            year = year + 1900;
            month = double(obj.get_scalar(9 * [1 1], 'uint8', filt) ) + 1;
            day = double(obj.get_scalar(10 * [1 1], 'uint8', filt) );
            hour = double(obj.get_scalar(11 * [1 1], 'uint8', filt) );
            minute = double(obj.get_scalar(12 * [1 1], 'uint8', filt) );
            second = double(obj.get_scalar(13 * [1 1], 'uint8', filt) );
            msec = double(obj.get_scalar(14 * [1 1], 'uint16', filt) );
            msec = msec / 10;
            val = datetime(year, month, day, hour, minute, second, msec,...
                'TimeZone',obj.timezone);
        end

        function val = get_distmidfirstcell(obj, filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            val = obj.get_cellsize(filt) + obj.get_blanking(filt);
        end
        function val = get_temperature(obj, filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            val = double(obj.get_scalar([18 18], 'int16', filt) );
            val = val / 100;
            val(~obj.burst_has_data(nortek.BurstBit.Temperature)) = nan;
        end
        function val = get_salinity(obj)
            val = nan(1, obj.nensembles);
        end
        function val = get_pressure(obj, filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            val = double(obj.get_scalar([20 20], 'uint32', filt) );
            val = val * 10;
            val(~obj.burst_has_data(nortek.BurstBit.Pressure)) = nan;
        end
        function val = get_echo(obj, filt)
            if nargin < 2
                filt = obj.get_data_filt;
            end
            if ~ obj.burst_has_data(nortek.BurstBit.Amplitude, filt)
                error('Amplitude data not available in burst record')
            end
            offs = obj.burst_data_offset(nortek.BurstBit.Amplitude, filt);
            val = double(obj.get_field(offs, 'uint8', filt))*.5;
        end
        function val = get_backscatter(obj)
        end
        function val = get_transducer(obj)
        end

        %%% METHODS AIDING READING OF DATA
        function dat = read_file(obj, file_name)
            if ischar(file_name)
                cur_files = dir(file_name);
                if isempty(cur_files)
                    error(['Could not find ' ,file_name])
                end
            elseif iscellstr(file_name)
                for cf = numel(file_name):-1:1
                    tmp_file = dir(file_name{cf});
                    if isempty(tmp_file)
                        error(['Could not find ', file_name{cf}])
                    end
                    cur_files(cf) = tmp_file;
                end
            end
            cur_files([cur_files.isdir]) = [];
            bytes = [cur_files.bytes];
            tot_bytes = sum(bytes);
            dat = zeros(tot_bytes, 1, 'uint8');
            cur_pos = 1;
            for cf = 1 : numel(cur_files)
                obj.files(cf) = string( ...
                    fullfile(cur_files(cf).folder, cur_files(cf).name) ...
                    );
                fid = fopen(obj.files(cf), 'r');
                end_pos = cur_pos + bytes(cf) - 1;
                [dat(cur_pos : end_pos), count] = fread(fid, '*uint8');
                fclose(fid);
                cur_pos = cur_pos + count;
            end
            if cur_pos < tot_bytes
                dat(cur_pos : end) = [];
            end
        end
        function offs = burst_data_offset(obj, bit, filt)
            if nargin < 3
                filt = obj.get_data_filt;
            end
            offs = zeros(1, obj.get_nensembles(filt));
            bv = obj.get_burst_version(filt);
            offs(bv == 2) = offs(bv == 2) + 68;
            offs(bv == 3 | bv == 1) = offs(bv == 3 | bv == 1) + 76;
            nc = obj.get_ncells(filt);
            nc_sounder = obj.get_sounder_ncells(filt);
            nfields = obj.get_nbeams(filt) .* nc;
            check_bit = @(x) bit > x & obj.burst_has_data(x, filt);
            offs = offs + check_bit(nortek.BurstBit.Velocity) .* nfields * 2;
            offs = offs + check_bit(nortek.BurstBit.Amplitude) .* nfields; 
            offs = offs + check_bit(nortek.BurstBit.Correlation) .* nfields;           
            offs = offs + check_bit(nortek.BurstBit.Altimeter) * 8;
            % for AST and Altimeter raw data, bit order does not match
            % order in data, so check are a bit more complicated.
            add_ast = (bit > nortek.BurstBit.AST |...
                bit == nortek.BurstBit.AltimeterRaw) &... 
                obj.burst_has_data(nortek.BurstBit.AST, filt);
            offs = offs + add_ast * 20;
            add_alt_raw = bit > nortek.BurstBit.AST & ...
                obj.burst_has_data(nortek.BurstBit.AltimeterRaw, filt);
            if bit > nortek.BurstBit.AST
                offs = offs + add_alt_raw .*...
                    (4 + 2 * obj.get_number_altimeter_samples(filt));
            end
            offs = offs + check_bit(nortek.BurstBit.EchoSounder) .*...
                nc_sounder * 2;
            offs = offs + check_bit(nortek.BurstBit.AHRS) * 16 * 4;
            offs = offs + check_bit(nortek.BurstBit.PercGood) .* nc;
        end

        function find_headers(obj)
            cur_pos = 0;
            c_ens = 1;
            est_size = floor(obj.nbytes / 1000);
            store_pos = nan(est_size, 1);
            store_head = zeros(est_size, 1, 'uint8');
            store_size = nan(est_size, 1);
            tot_bytes = obj.nbytes;
            while cur_pos < tot_bytes
                cur_pos = cur_pos + 1;
                f_sync = 1;
                if ~obj.raw(cur_pos) == nortek.DataHeader.AD2CP
                    f_sync = find(obj.raw(cur_pos : end) == ...
                        nortek.DataHeader.AD2CP, 1, 'first');
                end
                if isempty(f_sync)
                    break
                end
                f_sync = f_sync + cur_pos - 1;
                if f_sync + 1 > tot_bytes
                    break
                end
                head_size = double(obj.raw(f_sync + 1));
                if f_sync + head_size > tot_bytes
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
                size_data = double(typecast(obj.raw(...
                    f_sync + 4 : f_sync + 5), 'uint16'));
                if f_sync + head_size + size_data - 1 > tot_bytes
                    continue
                end
                csum_calc = obj.compute_checksum(f_sync + head_size,...
                    size_data);
                csum_read = typecast(obj.raw(f_sync + 6: f_sync + 7),...
                    'uint16');
                if csum_calc~=csum_read
                    continue
                end
                if c_ens > est_size
                    est_size = est_size*2;
                    store_pos(est_size) = nan;
                    store_head(est_size) = 0x00;
                    store_size(est_size) = nan;
                end
                store_pos(c_ens) = f_sync + head_size;
                store_head(c_ens) = obj.raw(f_sync + 2);
                store_size(c_ens) = size_data;
                c_ens = c_ens + 1;
                cur_pos = f_sync + head_size + size_data - 1;
            end
            store_pos(c_ens:end) = [];
            store_head(c_ens:end) = [];
            store_size(c_ens:end) = [];
            obj.data_position = store_pos;
            obj.data_header = store_head;
            obj.data_size =store_size;
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
        function filt = get_sounder_filt(obj, conf_id, sounder_id)
            if nargin < 2 || ismepty(conf_id)
                conf_id = obj.configuration_selector;
            end
            filt = obj.get_data_filt(nortek.DataHeader.EchoSounder, conf_id);
            if nargin < 3 || isempty(sounder_id)
                sounder_id = obj.echo_sounder_selector;
            end
            sound_filt = false(size(filt));
            sound_idx = reshape(obj.echo_sounder_index,[],1);
            for c_id = 1:numel(sounder_id)
                sound_filt = sound_filt | sound_idx == sounder_id(c_id);
            end
            filt = filt & sound_filt;
        end
        function filt = get_data_filt(obj, data_id, conf_id)
            if nargin < 2 || isempty(data_id)
                data_id = obj.data_selector;
            end
            filt = false(numel(obj.data_header),1);
            for c_id = 1:numel(data_id)
                filt = filt | obj.data_header == data_id(c_id);
            end
            if nargin < 3 || isempty(conf_id)
                conf_id = obj.configuration_selector;
            end
            filt_conf = false(numel(obj.data_header),1);
            act_conf = obj.active_configuration;
            for c_id = 1:numel(conf_id)
                c_filt = act_conf == conf_id(c_id);
                filt_conf = filt_conf | c_filt;
            end
            filt = filt & filt_conf;
        end
        function val = get_scalar(obj, offset, data_type, filt)
            if nargin < 4
                filt = obj.get_data_filt;
            end
            of = obj.data_position(filt);
            of = reshape(of, 1, []);
            ver = obj.get_burst_version(filt);
            if numel(offset) == 2
                ver3 = ver == 3 | ver == 1; % 1 is written in BT data
                ver2 = ver == 2;
                of(ver2) = of(ver2) + offset(1);
                of(ver3) = of(ver3) + offset(2);
            else
                if any(ver == 2)
                    error('Version 2 not supported by current data type')
                end
                of = of + offset(1);
            end
            of = (0 : obj.sizeof(data_type) -1 )' + of;
            of = of(:);
            val = typecast(obj.raw(of), data_type);
        end
        function val = get_field(obj, offset, data_type, filt)
            if nargin < 4
                filt = obj.get_data_filt;
            end
            of = obj.data_position(filt);
            of = reshape(of, 1, []) + offset;
            siz = obj.sizeof(data_type);
            nc = obj.get_ncells(filt);
            nb = obj.get_nbeams(filt);
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


    end
    properties (Access = protected)
        raw (1, :) uint8
        data_position(:, 1) double
        data_header(:,1) nortek.DataHeader
        data_size(:,1) double
    end
    properties (SetAccess = protected)
        files (:, 1) string
    end
    methods(Static)
        function val = get_bit(dat, start_idx, end_idx)
            % 0 based indices!!!
            dat_type = class(dat);
            if nargin < 3
                end_idx = start_idx;
            end
            start_idx = double(start_idx);
            end_idx = double(end_idx);
            val = bitand(bitshift(dat, - start_idx,dat_type), ...
                cast(2^(end_idx-start_idx+1),dat_type)-1, ...
                dat_type);
        end
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