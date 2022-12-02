classdef ADCP < ADCP
    % Wrapper class for adcp structures
    %
    %   obj=ADCP() Constructs default object
    %
    %   obj=ADCP(...) You can pass different objects to the class on
    %   initialization. Depending on the class of the object it will be
    %   assigned to a property:
    %   - Filter objects are appended to the filters property
    %   - acoustics.Water is assigned to the water property
    %   - acoustics.PistonTransducer is assigned to the transducer property
    %   - ADCPHorizontalPosition to horizontal_position_provider
    %   - ADCPVerticalPosition to vertical_position_provider
    %   - struct objects are assigned to the raw property 
    %
    %   ADCP properties:
    %   raw - adcp structure read by readADCP.m
    %   filters - filters for profiled data
    %   timezone - timezone of the data
    %   type - type of ADCP being used
    %   transformation_matrix_source - source of instrument matrix
    %   vertical_position_provider - Class providing vertical positions
    %   horizontal_position_provider- Class providing horizontal positions
    %   heading_provider - Class providing headings
    %
    %   ADCP read-only properties:
    %   *Ambient properties*
    %   temperature - ambient temperature (degrees Celsius)
    %   salinity - salinity (psu)
    %   pressure - pressure (Pa)
    %   water - acoustics.Water object describing water characteristics
    %
    %   *Instrument characteristics*
    %   transducer - object describing transducer characteristics
    %   is_workhorse - return whether ADCP is a Workhorse ADCP
    %
    %   *Sizing*
    %   fileid - ID of file ensemble was read from
    %   time - time ensembles were measured
    %   nensembles - number of ensembles
    %   nbeams - number of beams
    %   ncells - number of cells
    %
    %   *Coordinate systems properties*
    %   coordinate_system - coordinate system used
    %   beam_angle - angle of acoustic beam with vertical
    %   convexity - convexity of the instrument
    %   three_beam_solutions_used - whether three beams solutions were used
    %   tilts_used_in_transform - whether tilts were used in transformations
    %   bin_mapping_used - wheteher bin mapping was used during measurements
    %
    %   *Orientation*
    %   pitch - pitch rotation angle
    %   roll - roll rotation
    %   heading - heading rotation
    %   headalign - heading alignment
    %   is_upward - whether instrument was pointing upward
    %
    %   *Data positioning*
    %   cellsize - vertical size of depth cell
    %   lengthxmitpulse - length of transmitted pulse
    %   blanking - blanking distance
    %   distmidfirstcell - vertical distance to the first measured depth cell
    %   depth_cell_slant_range - slant range to depth cell
    %   horizontal_position - horizontal position of the ADCP
    %   vertical_position - vertical position of the ADCP
    %
    %   *Backscatter*
    %   bandwidth - bandwidth used for measurements (0=wide, 1=narrow)
    %   current - transmit current of transducer (A)
    %   current_factor - factor for current computation from ADC channel
    %   voltage - transmit voltage of transducer (V)
    %   voltage_factor - factor for voltage computation from ADC channel
    %   power - transmit power of transducer (W)
    %   attitude_temperature - temperature of transducer (Celsius)
    %   intensity_scale - intensity scale factor (dB/m)
    %   echo - Raw echo intennstiy (dB)
    %   backscatter_constant - Instrument constant (dB)
    %   backscatter - Volumne backscatter strength (dB)
    %
    %   ADCP methods
    %   bad - filter in use
    %   xform - coordinate transformation matrices
    %   depth_cell_offset - xyz offset from transducer to depth cell
    %   velocity - get velocity profile data
    %   plot_filters - plot active profile data filters
    %   plot_orientations - plot orientations of ADCP
    %   plot_velocity - plot velocity profiling data
    %
    %   ADCP static methods:
    %   beam_2_instrument - beam to instrument transformation matrices
    %   head_tilt - tilt matrices
    %   instrument_2_beam - instrument to beam transformation matrices
    %   int16_to_double - transform int16 to double values
    %   invert_xform - invert transformation matrices
    %
    % see also:VMADCP
    properties(SetObservable, AbortSet)
        % ADCP/raw
        %
        %   raw adcp structure as read with readADCP.m. Setting this
        %   property will reset the transducer and water properties.
        %
        % see also: ADCP, transducer, water
        raw(1,1) struct;
        
        % ADCP/temperature_offset
        %
        %   Specifies the temperature offset for Workhorse ADCPs for ADC
        %   channels. Default is -0.35. PS0 result
        %
        % see also: ADCP
        temperature_offset(1,1) double = -0.35
        
        


        on_battery(1,1) logical = false;
    end
    properties(Dependent, SetAccess=private)     
        % ADCP/type
        %
        %   Get the type of ADCP that is being used.
        %
        % see also: ADCP, rdi.ADCP_Type
        type rdi.ADCP_Type 


        % ADCP/fileid read only property
        %
        % fileid contains a fileid for each ensemble corresponding to the
        % file from which it was read
        %
        % see also: ADCP
        fileid
        
        % ADCP/convexity read only property
        %
        % convexity of the instrument. 1 for convex, -1 for concave
        %
        % see also: ADCP
        convexity
        
        % ADCP/is_upward read only property
        %
        % logical value indicating whether instrument is pointing upward.
        %
        % see also: ADCP
        is_upward
        
        % ADCP/tilts_used_in_transform read only property
        %
        % logical value indicating whether tilts were used in
        % transformation
        %
        % see also: ADCP
        tilts_used_in_transform
        
        % ADCP/bin_mapping_used read only property
        %
        % logical value indicating whether bin mapping was used.
        %
        % see also: ADCP
        bin_mapping_used
        
        % ADCP/three_beam_solutions_used  read only property
        %
        % logical value indicating whether three beam solutions were used
        %
        % see also: ADCP
        three_beam_solutions_used
        
        % ADCP/blanking  read only property
        %
        % blanking distance (m), i.e. unmeasured distance near transducer
        % due to ringing.
        %
        % see also: ADCP
        blanking
        
        % ADCP/lengthxmitpulse  read only property
        %
        % length of transmitted pulse (m), length of transmitted pulse
        %
        % see also: ADCP
        lengthxmitpulse
        
        % ADCP/headalign read only property
        %
        % heading alignment, i.e. offset between ship heading and ADCP
        % heading in degrees
        %
        % see also: ADCP
        headalign      

        % ADCP/bandwidth read only property
        %
        % bandwidth used (0=wide, 1=narrow)
        %
        % see also: ADCP
        bandwidth

        % ADCP/current read only property
        %
        % current on transducer (A)
        %
        % see also: ADCP
        current
        
        % ADCP/current_factor read only property
        %
        % factor to compute current from ADCP's ADC channels
        %
        % see also: ADCP
        current_factor
        
        % ADCP/voltage read only property
        %
        % voltage on transducer (V)
        %
        % see also: ADCP
        voltage
        
        % ADCP/voltage_factor read only property
        %
        % factor to convert ADC channel value to voltage
        %
        % see also: ADCP
        voltage_factor

        voltage_offset

        frequency
        
        % ADCP/power read only property
        %
        % Power on transducer (W)
        %
        % see also: ADCP
        power
        
        % ADCP/attitude_temperature read only property
        %
        % Attitude temperature of transducer (Celsius)
        %
        % see also: ADCP
        attitude_temperature
        
        % ADCP/intensity_scale
        %
        % Echo intensity scaling (dB/count)
        %
        % see also: ADCP
        intensity_scale
        
        % ADCP/backscatter_constant
        %
        % Instrument specific backscatter constant (dB)
        %
        % see also: ADCP
        backscatter_constant

        typical_pdbw

        % ADCP/noise_level
        %
        %   Specifies the background noise received by the instrument in
        %   counts. Default is -40. PT3 command
        %
        % see also: ADCP
        noise_level
    end

    properties(Constant)
        % frequency for Sentinel ADCPs
        sentinel_freq = [38.4, 76.8, 153.6, 307.2,...
            491.52, 983.04, 2457.6]*1e3

        % frequency for Workhorse ADCPs
        workhorse_freq = [38.4, 76.8, 153.6, 307.2,...
            614.4, 1228.8, 2457.6]*1e3

        % voltage scaling for workhorse ADCPs (from workhorse command and
        % output data format manual, Oct 2022, under ADC channels)
        wh_v_scaling = [2092719, 592157, 592157, 380667, 253765, 253765]

        % current scaling for workhorse ADCPs (from workhorse command and
        % output data format manual, Oct 2022, under ADC channels)
        wh_i_scaling = [43838, 11451, 11451, 11451, 11451, 11451, 11451]
    end

    methods
        %%% Constructor method %%%
        function obj=ADCP(varargin)
            obj = obj@ADCP(varargin{:});
            obj.instrument_matrix_provider = [
                rdi.InstrumentMatrixFromCalibration; 
                rdi.InstrumentMatrixFromBAngle
                ];
            obj.heading_provider = [
                rdi.HeadingProviderTFiles; 
                rdi.HeadingInternal
                ];
            for ca=1:nargin
                if isstruct(varargin{ca})
                    obj.raw=varargin{1};
                end
            end
        end
        
        %%% Set and Get methods %%%
        function file_id=get.fileid(obj)
            file_id=obj.raw.FileNumber;
        end
        function ha=get.headalign(obj)
            ha=double(obj.raw.headalign)/100;
        end
        function conv=get.convexity(obj)
            conv=reshape(bin2dec(obj.raw.sysconf(4,:)'),1,[]);
            conv(conv==0)=-1;
        end
        function up=get.is_upward(obj)
            up=reshape(bin2dec(obj.raw.sysconf(8,:)')==1,1,[]);
        end

        function tf=get.tilts_used_in_transform(obj)
            tf=reshape(bin2dec(obj.raw.corinfo(3,:)'),1,[])==1;
        end
        function tf=get.bin_mapping_used(obj)
            tf=reshape(bin2dec(obj.raw.corinfo(1,:)'),1,[])==1;
        end
        function tf=get.three_beam_solutions_used(obj)
            tf=reshape(bin2dec(obj.raw.corinfo(2,:)'),1,[])==1;
        end
        function blank=get.blanking(obj)
            blank=double(obj.raw.blnk)/100;
        end
        function l=get.lengthxmitpulse(obj)
            if isfield(obj.raw,'sp_transmit_length') % streampro leader support
                l=double(obj.raw.sp_transmit_length)/1000;
            else
                l=double(obj.raw.lngthtranspulse)/100;
            end
        end
        function val=get.bandwidth(obj)
            val=double(obj.raw.bandwidth);
        end
        function val = get.frequency(obj)
            val = nan(1,obj.nensembles);
            val(obj.type == rdi.ADCP_Type.STREAMPRO_31) = 2e6;
            val(obj.type == rdi.ADCP_Type.PINNACLE_61) = 44e3;
            id = bin2dec(fliplr(obj.raw.sysconf(1:3,:)'))'+2;
            is_sent = obj.type == rdi.ADCP_Type.SENTINELV_47 |...
                obj.type == rdi.ADCP_Type.SENTINELV_66;
            val(is_sent) = obj.sentinel_freq(id(is_sent));
            val(~is_sent) = obj.workhorse_freq(id(~is_sent));
        end
        function val = get.voltage_offset(obj)
            val = zeros(size(1,obj.nensembles));
            val(obj.type == rdi.ADCP_Type.PINNACLE_61)=24;
        end

        function val=get.voltage_factor(obj) % From workhorse operation manual and PDDecoder software
            val = 1e5*ones(size(obj.type));
            val(obj.type == rdi.ADCP_Type.TASMAN_74 |...
                obj.type == rdi.ADCP_Type.PATHFINDER_67) = 1e6;
            isw = obj.is_workhorse;
            freq = obj.frequency;
            [~, freq_idx] = ismember(freq, obj.workhorse_freq);
            fgood = freq_idx ~= 0;
            val(isw & fgood) = obj.wh_v_scaling(freq_idx(isw & fgood));
        end
        function val=get.current(obj)
            val=reshape(double(obj.raw.ADC(:,1)),1,[]).*obj.current_factor/1e6;
        end
        function val=get.voltage(obj)
            val=reshape(double(obj.raw.ADC(:,2)),1,[]).*...
                obj.voltage_factor/1e6+obj.voltage_offset;
        end
        function val=get.power(obj)
            %%% Has voltage and current
            has_volt = any(obj.raw.ADC(:,2)~=0);
            has_curr = any(obj.raw.ADC(:,1)~=0);
            val = nan(1,obj.nensembles);

            has_v_and_c = has_volt & has_curr;
            val(has_v_and_c) = obj.voltage(has_v_and_c) .* ...
                obj.current(has_v_and_c);

            
            %%% Has voltage, but no current
            only_v  = has_volt & ~has_curr;
            [~,pdbw_bat] = obj.get_instrument_characteristics(true);
            [~,pdbw_ps] = obj.get_instrument_characteristics(false);
            
            % workhorses
            iswh = obj.is_workhorse;
            val(iswh & only_v) = 10 .^ ( (pdbw_bat(iswh & only_v) +...
                20 * log10(obj.voltage(iswh & only_v) / 32) ) / 10);
            
            % pioneer and explorer
            ispe = obj.type == rdi.ADCP_Type.PIONEER_73 |...
                obj.type == rdi.ADCP_Type.EXPLORER_34;
            val(ispe & only_v) = 10 .^ ( (pdbw_ps(ispe & only_v) +...
                20 * log10(obj.voltage(ispe & only_v) / 32) ) / 10);

            % riogrande, channelmaster, riverpro, riverray
            isrg = obj.type == rdi.ADCP_Type.RIO_GRANDE_10 |...
                obj.type == rdi.ADCP_Type.CHANNELMASTER_28 |...
                obj.type == rdi.ADCP_Type.RIVERPRO_56 |...
                obj.type == rdi.ADCP_Type.RIVERRAY_44;
            val(isrg & only_v) = 10 .^ ( (pdbw_ps(isrg & only_v) +...
                20 * log10(obj.voltage(isrg & only_v) / 12) ) / 10);

            % sentinel V
            isst = obj.type == rdi.ADCP_Type.SENTINELV_47 |...
                obj.type == rdi.ADCP_Type.SENTINELV_66;
            val(isst & only_v) = 10 .^ ( (pdbw_bat(isst & only_v) +...
                20 * log10(obj.voltage(isst & only_v) / 14) ) / 10);

            %%% No voltage and no current
            typ = ~has_curr & ~has_volt;
            warning(['Make sure you set on_battery property correctly',...
                'for proper power computation.'])
            val(typ) = obj.typical_pdbw(typ);
        end
        function val=get.attitude_temperature(obj)
            if ~obj.is_workhorse
                warning('Assuming ADCP is a workhorse')
            end
            DC_COEF = 9.82697464e1;                                                    % Temperature coefficients
            FIRST_COEF = -5.86074151382e-3;
            SECOND_COEF = 1.60433886495e-7;
            THIRD_COEF = -2.32924716883e-12;
            t_cnts = reshape(double(obj.raw.ADC(:,6))*256,1,[]);                 % Temperature Counts (ADC value)
            val = obj.temperature_offset + ((THIRD_COEF.*t_cnts + SECOND_COEF).*t_cnts +...
                FIRST_COEF).*t_cnts + DC_COEF;                                         % real-time temperature of the transducer (C)
        end
        function val=get.intensity_scale(obj)
            val=127.3./(obj.attitude_temperature+273); % From WinRiver manual
            val(obj.type == rdi.ADCP_Type.RIVERRAY_44 |...
                obj.type == rdi.ADCP_Type.RIVERPRO_56 |...
                obj.type == rdi.ADCP_Type.EXPLORER_34 |...
                obj.type == rdi.ADCP_Type.PIONEER_73) = 0.6;
            val(obj.type == rdi.ADCP_Type.SENTINELV_47 |...
                obj.type == rdi.ADCP_Type.SENTINELV_66) = 0.5;
        end
        function val=get.noise_level(obj)
            val=ones(1,nensembles)*40;
            val(obj.type == rdi.ADCP_Type.RIVERRAY_44 |...
                obj.type == rdi.ADCP_Type.RIVERPRO_56 |...
                obj.type == rdi.ADCP_Type.EXPLORER_34 |...
                obj.type == rdi.ADCP_Type.PIONEER_73) = 50;
            val(obj.type == rdi.ADCP_Type.SENTINELV_47 |...
                obj.type == rdi.ADCP_Type.SENTINELV_66) = 40;
        end
        function val = get.typical_pdbw(obj)
            [~, val] = obj.get_instrument_characteristics;
        end

        function val=get.backscatter_constant(obj)
            val = obj.get_instrument_characteristics;
        end

        function val=get.type(obj)
            val = rdi.ADCP_Type(obj.raw.firmver);
        end
        %%% Ordinary methods
        
        function val=get.current_factor(obj) % From workhorse operation manual
            val = 1e5*ones(1,obj.nensembles);
            isw = obj.is_workhorse;
            freq = obj.frequency;
            [~, freq_idx] = ismember(freq, obj.workhorse_freq);
            fgood = freq_idx ~= 0;
            val(isw & fgood) = obj.wh_i_scaling(freq_idx(isw & fgood));
        end
        
        function tf=is_workhorse(obj)
           tf=ismember(obj.type,[...
               rdi.ADCP_Type.WORKHORSE_8;...
               rdi.ADCP_Type.WORKHORSE_16;...
               rdi.ADCP_Type.WORKHORSE_43;...
               rdi.ADCP_Type.WORKHORSE_50;...
               rdi.ADCP_Type.WORKHORSE_51;...
               rdi.ADCP_Type.WORKHORSE_52;...
               rdi.ADCP_Type.WORKHORSE_77;...
               rdi.ADCP_Type.WORKHORSE_78...                                                                      
               ]);
        end
        function vel=velocity(obj,dst,filter)
            vel=obj.int16_to_double(obj.raw.VEL)/1000;
            B=CoordinateSystem.Beam;
            src=obj.coordinate_system;
            if nargin > 1 && ~all(dst == src)
                if any(dst==B & ~(src==B) & obj.bin_mapping_used)
                    warning('Bin mapping was used: backward transformation to beam coordinates might be incorrect')
                end
                tm=obj.xform(dst);
                vel=helpers.matmult(tm, vel);
            end
            if nargin > 2
                bad=obj.bad(filter);
            else
                bad=obj.bad();
            end
            vel(bad)=nan;
        end
        function tm=xform(obj,dst, src,varargin)
            if nargin < 3
                src=obj.coordinate_system;
            end
            P=inputParser;
            P.addParameter('Geometry', false,@(x) isscalar(x) && islogical(x));
            P.parse(varargin{:});
            
            tm = repmat(shiftdim(eye(4),-2),1,obj.nensembles);
            I=CoordinateSystem.Instrument;
            S=CoordinateSystem.Ship;
            E=CoordinateSystem.Earth;
            B=CoordinateSystem.Beam;
            exp_cfilt=true(1,obj.nensembles); % dummy filter to expand scalar input
            croll=obj.roll;
            croll(obj.is_upward)=croll(obj.is_upward)+180;
            ha=obj.headalign;
            cpitch=obj.pitch;
            head=obj.heading;
            
            % FORWARD
            % from lower than instrument to instrument
            cfilt = exp_cfilt & dst >= I & src < I;
            tmptm=obj.instrument_matrix_provider.b2i_matrix(obj);
            tm(:,cfilt,:,:)=tmptm(:,cfilt,:,:);
            
            % from lower than ship to ship
            cfilt = exp_cfilt & dst >= S & src < S;
            tm(1,cfilt,:,:)=helpers.matmult(...
                obj.head_tilt(ha(cfilt),cpitch(cfilt), croll(cfilt)),...
                tm(1,cfilt,:,:));
            
            % from lower than earth to earth
            cfilt = exp_cfilt & dst >= E & src < E;
            tm(:,cfilt,:,:)=helpers.matmult(...
                obj.head_tilt(head(cfilt)),...
                tm(:,cfilt,:,:));
            
            % INVERSE
            % from higher than ship to ship
            cfilt = exp_cfilt & dst <= S & src > S;
            tm(:,cfilt,:,:)=permute(obj.head_tilt(head(cfilt)),[1,2,4,3]);
            
            % from higher than instrument to instrument
            cfilt = exp_cfilt & dst <= I & src > I;
            if any(strcmp(P.UsingDefaults,'Geometry'))
                cpitch(~obj.tilts_used_in_transform)=0; % Take into account whether tilts where used when transforming back
                croll(~obj.tilts_used_in_transform)=0; % Take into account whether tilts where used when transforming back
            elseif ~P.Results.Geometry
                cpitch(:)=0;
                croll(:)=0;
            end
                
            tm(1,cfilt,:,:)=helpers.matmult(...
                permute(obj.head_tilt(ha(cfilt),cpitch(cfilt), croll(cfilt)),[1,2,4,3]),...
                tm(1,cfilt,:,:));
            
            % from higher than beam to beam
            cfilt = exp_cfilt & dst == B & src > B;
            tmptm=obj.instrument_matrix_provider.i2b_matrix(obj);
            tm(:,cfilt,:,:)=helpers.matmult(...
                tmptm(:,cfilt,:,:),...
                tm(1,cfilt,:,:));
        end
    end
    methods(Access = protected)
        function cs=get_coordinate_system(obj)
            csnum=reshape(bin2dec(obj.raw.corinfo(4:5,:)'),1,[]);
            cs(1,obj.nensembles)=CoordinateSystem.Beam;
            cs(csnum==2)=CoordinateSystem.Instrument;
            cs(csnum==1)=CoordinateSystem.Ship;
            cs(csnum==3)=CoordinateSystem.Earth;
        end
        function ang=get_beam_angle(obj)
            ang=double(obj.raw.HADCPbeamangle);
            ang_sys=reshape(bin2dec(obj.raw.sysconf(9:10,:)'),1,[]);
            ang(ang==0 & ang_sys==0)=15;
            ang(ang==0 & ang_sys==2)=20;
            ang(ang==0 & ang_sys==3)=30;
            ang(ang==0)=nan;
        end
        function nens=get_nensembles(obj)
            nens=size(obj.fileid,2);
        end
        function nb=get_nbeams(obj)
            nb=double(obj.raw.usedbeams);
        end
        function c=get_cellsize(obj)
            if isfield(obj.raw,'sp_bin_space') %streampro leader support (space makes more sense than size, for the actual use. These seem always to match btw)
                c=double(obj.raw.sp_bin_space)/1000;
            else
                c=double(obj.raw.binsize)/100;
            end
        end
        function val=get_temperature(obj)
            val=double(obj.raw.temperature)/100;
        end
        function val=get_salinity(obj)
            val=double(obj.raw.salinity);
        end
        function val=get_pressure(obj)
            funderflow=obj.raw.pressure>3e9;
            val=double(obj.raw.pressure);
            val(funderflow)=val(funderflow)-double(intmax('uint32'));
            val=val*10;
        end
        function t=get_time(obj)
            t=reshape(datetime(obj.raw.timeV,'TimeZone',obj.timezone),1,[]);
        end
        function db1=get_distmidfirstcell(obj)
            if isfield(obj.raw,'sp_mid_bin1') % streampro leader support (more accurate)
                db1=double(obj.raw.sp_mid_bin1)/1000;
            else
                db1=double(obj.raw.distmidbin1)/100;
            end
        end

        function n=get_ncells(obj)
            n=double(obj.raw.nbins);
        end
        function val=get_echo(obj)
            val=double(obj.raw.ECHO).*obj.intensity_scale;
        end

        function [C, Pdbw, rayl] = get_instrument_characteristics(obj, bat)
            % characteristics from FS031
            tp = obj.type;
            freq = obj.frequency;
            isw = obj.is_workhorse;
            wh_freq = obj.workhorse_freq;
            st_freq = obj.sentinel_freq;
            bw = obj.bandwidth;
            if nargin < 2
                bat = obj.on_battery;
            end
            C = nan(1,obj.nensembles);
            Pdbw = nan(1,obj.nensembles);
            rayl = nan(1,obj.nensembles);

            iscm = tp == rdi.ADCP_Type.CHANNELMASTER_28;
            C(iscm & freq == wh_freq(4) & bw == 0) = -143.4; % 300 high
            C(iscm & freq == wh_freq(4) & bw == 1) = -152.26; % 300 low
            C(iscm & freq == wh_freq(5) & bw == 0) = -139.08; % 600 high
            C(iscm & freq == wh_freq(5) & bw == 1) = -147.28; % 600 low
            C(iscm & freq == wh_freq(6) & bw == 0) = -127.13; % 1200 high
            C(iscm & freq == wh_freq(6) & bw == 1) = -137.17; % 1200 low

            Pdbw(iscm & freq == wh_freq(4)) = 15.1; % 300 high
            Pdbw(iscm & freq == wh_freq(5)) = 12.0; % 600 high
            Pdbw(iscm & freq == wh_freq(6)) = 9.0; % 1200 high

            rayl(iscm & freq == wh_freq(4)) = 2.69; % 300 high
            rayl(iscm & freq == wh_freq(5)) = 2.96; % 600 high
            rayl(iscm & freq == wh_freq(6)) = 1.71; % 1200 high


            isex = tp == rdi.ADCP_Type.EXPLORER_34;
            if any(isex)
                obj.warning(...
                    'Assuming Explorer ADCP with piston transducer')
            end
            C(isex & bw == 0) = -132.73; % high
            C(isex & bw == 1) = -140.95; % low
            Pdbw(isex) = 3;
            rayl(isex) = 1.35;
            
            % long ranger (only workhorse with 75 KHz)
            islr = isw & freq == wh_freq(2);
            C(islr & bw == 0) = -161.19;
            C(islr & bw == 1) = -166.94;
            Pdbw(islr & bat) = 23.8;
            Pdbw(islr & ~bat) = 27.3;
            rayl(islr) = 1.26;


            isos = tp == rdi.ADCP_Type.OCEAN_SURVEYOR_14 |...
                tp == rdi.ADCP_Type.OCEAN_SURVEYOR_23;
            C(isos & freq == wh_freq(1)) = -172.19; % 38
            C(isos & freq == wh_freq(2)) = -164.26; % 75
            C(isos & freq == wh_freq(3)) = -156.01; % 150
            Pdbw(isos & freq == wh_freq(1)) = 24.0; % 38
            Pdbw(isos & freq == wh_freq(2)) = 24.0; % 75
            Pdbw(isos & freq == wh_freq(3)) = 21.0; % 150
            rayl(isos & freq == wh_freq(1)) = 8.19; % 38
            rayl(isos & freq == wh_freq(2)) = 3.24; % 75
            rayl(isos & freq == wh_freq(3)) = 1.62; % 150

            ispi = tp == rdi.ADCP_Type.PIONEER_73;
            C(ispi & freq == wh_freq(4)) = -151.30; % 300
            C(ispi & freq == wh_freq(5)) = -145.25; % 600
            Pdbw(ispi & freq == wh_freq(4)) = 14; % 300
            Pdbw(ispi & freq == wh_freq(5)) = 9; % 600
            rayl(ispi & freq == wh_freq(4)) = 1.67; % 300
            rayl(ispi & freq == wh_freq(5)) = 1.31; % 600
            
            % quarter master (only workhorse on 150 KHz)
            isqm = isw & freq == wh_freq(3);
            C(isqm & bw == 0) = -153.75;
            C(isqm & bw == 1) = -161.01;
            Pdbw(isqm & bat) = 15.1;
            Pdbw(isqm & ~bat) = 18.6;
            rayl(isqm) = 1.68;
            

            isrg = tp == rdi.ADCP_Type.RIO_GRANDE_10;
            C(isrg & freq == wh_freq(5) & bw == 0) = -139.09; % 600 high
            C(isrg & freq == wh_freq(5) & bw == 1) = -149.14; % 600 low
            C(isrg & freq == wh_freq(6) & bw == 0) = -129.44; % 1200 high
            C(isrg & freq == wh_freq(6) & bw == 1) = -139.57; % 1200 low
            Pdbw(isrg & freq == wh_freq(5)) = 9; % 600
            Pdbw(isrg & freq == wh_freq(6)) = 4.8; % 1200
            rayl(isrg & freq == wh_freq(5)) = 1.75; % 600
            rayl(isrg & freq == wh_freq(6)) = 1.71; % 1200

            %Riverpro (5b) or RioPro (4b)
            isrp = tp == rdi.ADCP_Type.RIVERPRO_56;
            beamconf = bin2dec(obj.raw.sysconf(13:16,:)')';
            is5b = beamconf == 10 | beamconf == 15;
            % RioPro
            C(isrp & ~is5b & freq == wh_freq(6) & bw == 0) = -131.36;
            C(isrp & ~is5b & freq == wh_freq(6) & bw == 1) = -141.08;
            Pdbw(isrp & ~is5b & freq == wh_freq(6)) = 7.8;
            rayl(isrp & ~is5b & freq == wh_freq(6)) = 1.71;
            
            % RiverPro
            C(isrp & is5b & freq == wh_freq(6) & bw == 0) = -128.09;
            C(isrp & is5b & freq == wh_freq(6) & bw == 1) = -137.81;
            Pdbw(isrp & is5b & freq == wh_freq(6)) = 7.8;
            rayl(isrp & is5b & freq == wh_freq(6)) = 0.81;

            isrr = tp == rdi.ADCP_Type.RIVERRAY_44;
            C(isrr) = -138.02;
            Pdbw(isrr) = 9;
            rayl(isrr) = 1.31;

            isst = tp == rdi.ADCP_Type.SENTINELV_47 |...
                tp == rdi.ADCP_Type.SENTINELV_66;
            % V100: 300 KHz, V50: 500 KHz, V100: 1000 KHz
            C(isst & freq == st_freq(4) & bw == 0) = -144.74; %V100 high
            C(isst & freq == st_freq(4) & bw == 1) = -151.24; %V100 low
            C(isst & freq == st_freq(5) & bw == 0) = -139.18; %V50 high
            C(isst & freq == st_freq(5) & bw == 1) = -145.73; %V50 low
            C(isst & freq == st_freq(6) & bw == 0) = -135.49; %V20 high
            C(isst & freq == st_freq(6) & bw == 1) = -143.32; %V20 low
            Pdbw(isst & freq == st_freq(4) & bat) = 14; %V100 bat
            Pdbw(isst & freq == st_freq(4) & ~bat) = 16.2; %V100 ps
            Pdbw(isst & freq == st_freq(5) & bat) = 10.8; %V50 bat
            Pdbw(isst & freq == st_freq(5) & ~bat) = 13; %V50 ps
            Pdbw(isst & freq == st_freq(6) & bat) = 9; %V20 bat
            Pdbw(isst & freq == st_freq(6) & ~bat) = 11.2; %V20 ps
            rayl(isst & freq == st_freq(4)) = 0.86; %V100 bat
            rayl(isst & freq == st_freq(4)) = 1.89; %V100 ps
            rayl(isst & freq == st_freq(5)) = 1.22; %V50 bat

            % workhorse 
            C(isw & freq == wh_freq(4) & bw == 0) = -140.78; %WH300 high
            C(isw & freq == wh_freq(4) & bw == 1) = -151.64; %WH300 low
            C(isw & freq == wh_freq(5) & bw == 0) = -139.09; %WH600 high
            C(isw & freq == wh_freq(5) & bw == 1) = -149.14; %WH600 low
            C(isw & freq == wh_freq(6) & bw == 0) = -129.44; %WH1200 high
            C(isw & freq == wh_freq(6) & bw == 1) = -139.57; %WH600 low
            Pdbw(isw & freq == wh_freq(4) & bat) = 14; %WH300 bat
            Pdbw(isw & freq == wh_freq(4) & ~bat) = 17.5; %WH300 ps
            Pdbw(isw & freq == wh_freq(5) & bat) = 9; %WH600 bat
            Pdbw(isw & freq == wh_freq(5) & ~bat) = 12.5; %WH600 ps
            Pdbw(isw & freq == wh_freq(6) & bat) = 4.8; %WH1200 bat
            Pdbw(isw & freq == wh_freq(6) & ~bat) = 8.3; %WH600 ps
            rayl(isw & freq == wh_freq(4)) = 0.87; %WH300
            rayl(isw & freq == wh_freq(5)) = 1.75; %WH600
            rayl(isw & freq == wh_freq(6)) = 1.71; %WH1200
        end

        function pt = get_transducer(obj, varargin)
        % Reset tranducer properties
        %
        %   reset_tranducer(obj) assigns the frequency and radius
        %   properties of the acoustics.Transducer object stored in the
        %   transducer property
        %
        % see also: ADCP, transducer
            
            pt=acoustics.PistonTransducer;
            %TODO: compute SentinelV and MonitorV radii back from
            %beam_width and frequency (see equation in Deines 1999)
            % Handle phased array ADCPs
            if obj.type==rdi.ADCP_Type.RIVERRAY_44
                if ~isa(pt,'acoustics.PhasedArrayTransducer')
                    pt=acoustics.PhasedArrayTransducer;
                end
                pt.frequency=614.4e3; % from RiverRay manual
                pt.radius=0.076/2; % from RiverRay manual
                return
            end
            
            % Piston transducer ADCPs
            if ~isa(pt,'acoustics.PistonTransducer')
                pt=acoustisc.PistonTransducer;
            end
            pt.water=obj.water;
            pt.depth=max(obj.pressure./9.81./pt.water.density,0);
            sysid=obj.raw.sysconf(1,1:3);
            switch sysid
                case '000' % Only Long Ranger ADCPs
                    pt.radius=0.203/2; % From operation Manual Long Ranger (drawing 6021, 6055)
                    pt.frequency=76.8e3; % From operation Manual Long Ranger
                case '100' % Only QuarterMaster ADCPs
                    pt.frequency=153.6e3; % From operation Manual Quarter Master
                    switch obj.type
                        case {rdi.ADCP_Type.QuarterMaster1500,rdi.ADCP_Type.QuarterMaster3000}
                            pt.radius=0.178/2; % From operation manual Quarter Master (drawing 6082, 6083)
                        case rdi.ADCP_Type.QuarterMaster1500ModBeams
                            pt.radius=0.1854/2; % From operation manual Quarter Master (drawing 1106)
                        case rdi.ADCP_Type.QuarterMaster6000
                            pt.radius=0.184/2; % From operation manual Quarter Master (drawing 1082)
                        otherwise
                            warning('Unknown radius for given ADCP type, assuming QM1500_Modular Beams')
                            pt.radius=0.1854/2; % From operation manual Quarter Master (drawing 1106)
                    end
                case '010'
                    pt.frequency=307.2e3; % Workhorse and SentinelV manual
                    switch obj.type
                        case {rdi.ADCP_Type.SentinelV,rdi.ADCP_Type.MonitorV}
                            pt.radius=nan;
                        case {rdi.ADCP_Type.Monitor,rdi.ADCP_Type.Sentinel}
                            pt.radius=0.0984/2; % Workhorse operation manual
                        case rdi.ADCP_Type.Mariner
                            pt.radius=0.0895/2; % Workhorse operation manual
                        otherwise
                            warning('Unknown frequency and radius for ADCP type, assuming Monitor or Sentinel')
                            pt.radius=0.0984/2; % Workhorse operation manual
                    end
                case '110'
                    switch obj.type
                        case {rdi.ADCP_Type.SentinelV}
                            pt.frequency=491.52e3; % Sentinel V operation manual
                            pt.radius=nan;
                        case {rdi.ADCP_Type.Monitor,rdi.ADCP_Type.Sentinel}
                            pt.frequency=614.4e3; % Workhorse operation manual
                            pt.radius=0.0984/2; % Workhorse operation manual
                        case {rdi.ADCP_Type.Mariner,rdi.ADCP_Type.RioGrande}
                            pt.frequency=614.4e3; % Workhorse operation manual
                            pt.radius=0.0895/2; % Workhorse operation manual
                        otherwise
                            warning('Unknown frequency and radius for ADCP type, assuming Monitor or Sentinel')
                            pt.frequency=614.4e3; % Workhorse operation manual
                            pt.radius=0.0984/2; % Workhorse operation manual
                    end
                case '001'
                    switch obj.type
                        case {rdi.ADCP_Type.SentinelV, rdi.ADCP_Type.MonitorV}
                            pt.frequency=983.04e3; % Sentinel V operation manual
                            pt.radius=nan;
                        case {rdi.ADCP_Type.Monitor , rdi.ADCP_Type.Sentinel, rdi.ADCP_Type.Mariner, rdi.ADCP_Type.RioGrande}
                            pt.frequency=1228.8e3; % Workhorse and RioGrande operation manual
                            pt.radius=0.0699/2; % Workhorse and RioGrande operation manual
                        otherwise
                            warning('Unknown frequency and radius for ADCP type, assuming Monitor, Mariner, Sentinel or Rio Grande')
                            pt.frequency=1228.8e3; % Workhorse and RioGrande operation manual
                            pt.radius=0.0699/2; % Workhorse and RioGrande operation manual
                    end
                case '101' % 2457.6e3 had this number, but don't know which ADCP has such specs?
                    pt.radius=nan;
                    switch obj.type
                        case rdi.ADCP_Type.StreamPro
                            pt.frequency=2e6; % From StreamPro Manual
                            pt.radius=nan;
                        otherwise
                            warning('Unknown frequency and radius for ADCP type, assuming StreamPro')
                            pt.frequency=2e6; % From StreamPro Manual
                            pt.radius=nan;
                    end
            end
        end
        function val=get_backscatter(obj)
            pt=obj.transducer;
            R=obj.depth_cell_slant_range+obj.cellsize/2/cosd(obj.beam_angle); % slant range to last quarter of cell
            two_alpha_R = 2.*pt.attenuation.*R;                    % compute 2alphaR
            LDBM=10*log10(obj.lengthxmitpulse);
            PDBW=10*log10(obj.power);                        
            val = obj.backscatter_constant + 10*log10((obj.attitude_temperature+273.16).*R.^2.*pt.near_field_correction(R).^2) - LDBM - PDBW + two_alpha_R + 10*log10(10.^((obj.echo-obj.noise_level)/10)-1); % equation according to fsa-031, correcting goustiaus and van haren equation      
            val(obj.bad)=nan;
        end
    
    end
    
    methods (Static)
        function tm=head_tilt(heading,pitch,roll)
            % Return the transformation matrix for the three tilts
            %
            %   tm=head_tilt(heading,pitch,roll) computes the rotation matrices
            %   given the heading, pitch and roll angles in degrees.
            %
            %   see also: ADCP, xform
            zr=zeros(size(heading));
            on=ones(size(heading));
            sh=sind(heading);
            ch=cosd(heading);
            if nargin < 2
                pitch = zr;
            end
            sp=sind(pitch);
            cp=cosd(pitch);
            if nargin < 3
                roll = zr;
            end
            sr=sind(roll);
            cr=cosd(roll);
            tm = cat(3,...
                cat(4,  ch.*cr+sh.*sp.*sr,  sh.*cp, ch.*sr-sh.*sp.*cr, zr),...
                cat(4, -sh.*cr+ch.*sp.*sr,  ch.*cp, -sh.*sr-ch.*sp.*cr, zr),...
                cat(4,            -cp.*sr,      sp,            cp.*cr, zr),...
                cat(4,                 zr,      zr,                zr, on));
        end
        function d=int16_to_double(i)
            % Transform 16 bit signed integers to double
            %
            %   d=int16_to_double(i) transforms the 16 bit integers in i to
            %   doubles in d. Integers equal to -32768 are assigned to be NaNs.
            fbad=i==-32768;
            d=double(i);
            d(fbad)=nan;
        end
        
    end
end