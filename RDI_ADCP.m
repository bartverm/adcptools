classdef RDI_ADCP < ADCP
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
        
        
        % ADCP/transformation_matrix_source
        %
        %   Specifies the sources for the transformation matrix as a
        %   InstrumentMatrixProvider. This is a vector which allows to
        %   define different sources. The first object capable of providing
        %   the matrix will be used.
        %
        % see also: ADCP, InstrumentMatrixProvider
        transformation_matrix_source (:,1) InstrumentMatrixProvider = [InstrumentMatrixFromCalibration; InstrumentMatrixFromBAngle];
        
        % ADCP/type
        %
        %   Specify the type of ADCP that is being used. Default is
        %   ADCP_Type.Monitor.
        %
        % see also: ADCP, ADCP_Type
        type(1,1) ADCP_Type = ADCP_Type.Unknown
        
        % ADCP/temperature_offset
        %
        %   Specifies the temperature offset for Workhorse ADCPs for ADC
        %   channels. Default is -0.35. PS0 result
        %
        % see also: ADCP
        temperature_offset(1,1) double = -0.35
        
        
        % ADCP/noise_level
        %
        %   Specifies the background noise received by the instrument in
        %   counts. Default is -40. PT3 command
        %
        % see also: ADCP
        noise_level(1,1) double = -40
    end
    properties(Dependent, SetAccess=private)
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
    end

    methods
        %%% Constructor method %%%
        function obj=RDI_ADCP(varargin)
            obj = obj@ADCP(varargin{:});
            obj.type=ADCP_Type.Unknown;
            for ca=1:nargin
                if isa(varargin{ca},'ADCP_Type')
                    obj.type=varargin{ca};
                elseif isstruct(varargin{ca})
                    obj.raw=varargin{1};
                end
            end
        end
        
        %%% Set and Get methods %%%
        function file_id=get.fileid(obj)
            file_id=obj.raw.FileNumber;
        end
        function ha=get.headalign(obj)
            ha=double(obj.raw.headalign(obj.fileid))/100;
        end
        function conv=get.convexity(obj)
            conv=reshape(bin2dec(obj.raw.sysconf(obj.fileid,4)),1,[]);
            conv(conv==0)=-1;
        end
        function up=get.is_upward(obj)
            up=reshape(bin2dec(obj.raw.sysconf(obj.fileid,8))==1,1,[]);
        end

        function tf=get.tilts_used_in_transform(obj)
            tf=reshape(bin2dec(obj.raw.corinfo(obj.fileid,3)),1,[])==1;
        end
        function tf=get.bin_mapping_used(obj)
            tf=reshape(bin2dec(obj.raw.corinfo(obj.fileid,1)),1,[])==1;
        end
        function tf=get.three_beam_solutions_used(obj)
            tf=reshape(bin2dec(obj.raw.corinfo(obj.fileid,2)),1,[])==1;
        end
        function blank=get.blanking(obj)
            blank=double(obj.raw.blnk(obj.fileid))/100;
        end
        function l=get.lengthxmitpulse(obj)
            if isfield(obj.raw,'sp_transmit_length') % streampro leader support
                l=double(obj.raw.sp_transmit_length)/1000;
            else
                l=double(obj.raw.lngthtranspulse(obj.fileid))/100;
            end
        end
        function val=get.bandwidth(obj)
            val=double(obj.raw.bandwidth(obj.fileid));
        end

        function val=get.voltage_factor(obj) % From workhorse operation manual
            if isfield(obj.raw,'sp_mid_bin1') % streampro header, use voltage factor from manual (0.1 to convert to voltage, but factor is divide by 1e6 in voltage calculation)
                val=1e5;
                return
            end
            if ~obj.is_workhorse
                warning('Assuming ADCP is a Workhorse')
            end
            switch obj.transducer.frequency
                case 76.8e3
                    val=2092719;
                case {153.6e3, 307.2e3}
                    val=592157;
                case 614.4e3
                    val=380667;
                case {1228.8e3, 2457.6e3}
                    val=253765;
                otherwise
                    warning('Unknown voltage factor')
                    val=nan;
            end
        end
        function val=get.current(obj)
            val=reshape(double(obj.raw.ADC(:,1))*obj.current_factor/1e6,1,[]);
        end
        function val=get.voltage(obj)
            val=reshape(double(obj.raw.ADC(:,2))*obj.voltage_factor/1e6,1,[]);
        end
        function val=get.power(obj)
            val=obj.voltage.*obj.current;
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
        end
        function val=get.backscatter_constant(obj)
            if ~is_workhorse(obj)
                warning('Assuming ADCP is a workhorse')
            end
            switch obj.transducer.frequency
                case 76.8e3
                    val=-159.1;
                case 307.2e3
                    switch obj.type
                        case ADCP_Type.Sentinel
                            val=-143.5;
                        case ADCP_Type.Monitor
                            val=-143;
                        otherwise
                            warning('Assuming ADCP is a Monitor')
                            val=-143;
                    end
                case 614.4e3
                    if obj.type~=ADCP_Type.RioGrande
                        warning('Assuming ADCP is a RioGrande')
                    end
                    val=-139.3;
                case 1228.8e3
                    val=-129.1;
                otherwise
                    warning('Unknown backscatter constant for current ADCP Type')
                    val=nan;
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
        
        %%% Ordinary methods
        
        function val=get.current_factor(obj) % From workhorse operation manual
            if ~obj.is_workhorse
                warning('Assuming ADCP is a Workhorse')
            end
            switch obj.transducer.frequency
                case 76.8e3
                    val=43838;
                case {153.6e3,307.2e3, 614.4e3, 1228.8e3, 2457.6e3}
                    val=11451;
                otherwise
                    warning('Do not know current factor for given ADCP type')
                    val=nan;
            end
        end
        
        function tf=is_workhorse(obj)
            if any(obj.type== [ADCP_Type.LongRanger1500, ADCP_Type.LongRanger3000, ADCP_Type.QuarterMaster1500,...
                    ADCP_Type.QuarterMaster1500ModBeams, ADCP_Type.QuarterMaster3000, ADCP_Type.QuarterMaster6000,...
                    ADCP_Type.Sentinel, ADCP_Type.Mariner, ADCP_Type.Monitor])
                tf=true;
            else
                tf=false;
            end
        end
        function vel=velocity(obj,dst,filter)
            % velocity profile data
            %
            %   vel=velocity(obj) returns the profiled velocity in m/s.
            %
            %   vel=velocity(obj,dst) specify destination coordinate system as
            %   CoordinateSystem object.
            %
            %   vel=velocity(obj,dst,filter) specify filter to be used insted
            %   of the ones specified in the current object.
            %
            %   see also: ADCP, CoordinateSystem
            vel=ADCP.int16_to_double(obj.raw.VEL)/1000;
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
            % Get transformation matrices for coordinate transformations
            %
            %   tm=xform(obj,dst) get the transformation matrices for the
            %   transformation from obj.coordinate_system to the given
            %   destination coordinate system. matrix will be 1 x nensembles x
            %   4 x 4
            %
            %   tm=xform(obj,dst,src) to specify a custom source coordinate
            %   system
            %
            %   tm=xform(...,'ParameterName', parameterValue) allows to
            %   specify the following options:
            %   'UseTilts' - if unset checks the ADCP data in inverse
            %   transformations to know whether to use the tilts. If set to
            %   true will always use the tilts in invers transformations,
            %   if set to false will never use tilts in inverse
            %   transformations.
            %
            %   see also: ADCP
            if nargin < 3
                src=obj.coordinate_system;
            end
            P=inputParser;
            P.addParameter('UseTilts', false,@(x) isscalar(x) && islogical(x));
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
            tmptm=obj.transformation_matrix_source.b2i_matrix(obj);
            tm(:,cfilt,:,:)=tmptm(:,cfilt,:,:);
            
            % from lower than ship to ship
            cfilt = exp_cfilt & dst >= S & src < S;
            tm(1,cfilt,:,:)=helpers.matmult(...
                ADCP.head_tilt(ha(cfilt),cpitch(cfilt), croll(cfilt)),...
                tm(1,cfilt,:,:));
            
            % from lower than earth to earth
            cfilt = exp_cfilt & dst >= E & src < E;
            tm(:,cfilt,:,:)=helpers.matmult(...
                ADCP.head_tilt(head(cfilt)),...
                tm(:,cfilt,:,:));
            
            % INVERSE
            % from higher than ship to ship
            cfilt = exp_cfilt & dst <= S & src > S;
            tm(:,cfilt,:,:)=permute(ADCP.head_tilt(head(cfilt)),[1,2,4,3]);
            
            % from higher than instrument to instrument
            cfilt = exp_cfilt & dst <= I & src > I;
            if any(strcmp(P.UsingDefaults,'UseTilts'))
                cpitch(~obj.tilts_used_in_transform)=0; % Take into account whether tilts where used when transforming back
                croll(~obj.tilts_used_in_transform)=0; % Take into account whether tilts where used when transforming back
            elseif ~P.Results.UseTilts
                cpitch(:)=0;
                croll(:)=0;
            end
                
            tm(1,cfilt,:,:)=helpers.matmult(...
                permute(ADCP.head_tilt(ha(cfilt),cpitch(cfilt), croll(cfilt)),[1,2,4,3]),...
                tm(1,cfilt,:,:));
            
            % from higher than beam to beam
            cfilt = exp_cfilt & dst == B & src > B;
            tmptm=obj.transformation_matrix_source.i2b_matrix(obj);
            tm(:,cfilt,:,:)=helpers.matmult(...
                tmptm(:,cfilt,:,:),...
                tm(1,cfilt,:,:));
        end
    end
    methods(Access = protected)
        function cs=get_coordinate_system(obj)
            csnum=reshape(bin2dec(obj.raw.corinfo(obj.fileid,4:5)),1,[]);
            cs(1,obj.nensembles)=CoordinateSystem.Beam;
            cs(csnum==2)=CoordinateSystem.Instrument;
            cs(csnum==1)=CoordinateSystem.Ship;
            cs(csnum==3)=CoordinateSystem.Earth;
        end
        function ang=get_beam_angle(obj)
            ang=double(obj.raw.HADCPbeamangle(obj.fileid));
            ang_sys=reshape(bin2dec(obj.raw.sysconf(obj.fileid,9:10)),1,[]);
            ang(ang==0 & ang_sys==0)=15;
            ang(ang==0 & ang_sys==2)=20;
            ang(ang==0 & ang_sys==3)=30;
            ang(ang==0)=nan;
        end
        function pitch=get_pitch(obj)
            pitch=double(obj.raw.pitch)/100;
            pitch=atand(tand(pitch).*cosd(obj.roll));
        end
        function roll=get_roll(obj)
            roll=double(obj.raw.roll)/100;
        end
        function head=get_heading(obj)
            head=obj.heading_provider.heading(obj);
        end
        function nens=get_nensembles(obj)
            nens=size(obj.fileid,2);
        end
        function nb=get_nbeams(obj)
            nb=double(obj.raw.usedbeams(obj.fileid));
        end
        function c=get_cellsize(obj)
            if isfield(obj.raw,'sp_bin_space') %streampro leader support (space makes more sense than size, for the actual use. These seem always to match btw)
                c=double(obj.raw.sp_bin_space)/1000;
            else
                c=double(obj.raw.binsize(obj.fileid))/100;
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
                db1=double(obj.raw.distmidbin1(obj.fileid))/100;
            end
        end

        function n=get_ncells(obj)
            n=double(obj.raw.nbins(obj.fileid));
        end
        function val=get_echo(obj)
            val=double(obj.raw.ECHO).*obj.intensity_scale;
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
            if obj.type==ADCP_Type.RiverRay
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
                        case {ADCP_Type.QuarterMaster1500,ADCP_Type.QuarterMaster3000}
                            pt.radius=0.178/2; % From operation manual Quarter Master (drawing 6082, 6083)
                        case ADCP_Type.QuarterMaster1500ModBeams
                            pt.radius=0.1854/2; % From operation manual Quarter Master (drawing 1106)
                        case ADCP_Type.QuarterMaster6000
                            pt.radius=0.184/2; % From operation manual Quarter Master (drawing 1082)
                        otherwise
                            warning('Unknown radius for given ADCP type, assuming QM1500_Modular Beams')
                            pt.radius=0.1854/2; % From operation manual Quarter Master (drawing 1106)
                    end
                case '010'
                    pt.frequency=307.2e3; % Workhorse and SentinelV manual
                    switch obj.type
                        case {ADCP_Type.SentinelV,ADCP_Type.MonitorV}
                            pt.radius=nan;
                        case {ADCP_Type.Monitor,ADCP_Type.Sentinel}
                            pt.radius=0.0984/2; % Workhorse operation manual
                        case ADCP_Type.Mariner
                            pt.radius=0.0895/2; % Workhorse operation manual
                        otherwise
                            warning('Unknown frequency and radius for ADCP type, assuming Monitor or Sentinel')
                            pt.radius=0.0984/2; % Workhorse operation manual
                    end
                case '110'
                    switch obj.type
                        case {ADCP_Type.SentinelV}
                            pt.frequency=491.52e3; % Sentinel V operation manual
                            pt.radius=nan;
                        case {ADCP_Type.Monitor,ADCP_Type.Sentinel}
                            pt.frequency=614.4e3; % Workhorse operation manual
                            pt.radius=0.0984/2; % Workhorse operation manual
                        case {ADCP_Type.Mariner,ADCP_Type.RioGrande}
                            pt.frequency=614.4e3; % Workhorse operation manual
                            pt.radius=0.0895/2; % Workhorse operation manual
                        otherwise
                            warning('Unknown frequency and radius for ADCP type, assuming Monitor or Sentinel')
                            pt.frequency=614.4e3; % Workhorse operation manual
                            pt.radius=0.0984/2; % Workhorse operation manual
                    end
                case '001'
                    switch obj.type
                        case {ADCP_Type.SentinelV, ADCP_Type.MonitorV}
                            pt.frequency=983.04e3; % Sentinel V operation manual
                            pt.radius=nan;
                        case {ADCP_Type.Monitor , ADCP_Type.Sentinel, ADCP_Type.Mariner, ADCP_Type.RioGrande}
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
                        case ADCP_Type.StreamPro
                            pt.frequency=2e6; % From StreamPro Manual
                            pt.radius=nan;
                        otherwise
                            warning('Unknown frequency and radius for ADCP type, assuming StreamPro')
                            pt.frequency=2e6; % From StreamPro Manual
                            pt.radius=nan;
                    end
            end
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
        function inv=invert_xform(tm)
            % Inverts transformation matrices
            %
            %   Inverts transformation matrices given in tm. Matrices are
            %   defined in the third and fourth dimension.
            %
            % see also: ADCP, xform
            inv=nan(size(tm));
            for ce=1:size(tm,2)
                inv(1,ce,:,:)=shiftdim(inv(squeeze(tm(1,ce,:,:))),-2);
            end
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