classdef PD0 < handle
% RDI PD0 binary format parser
% 
% PD0 Properties:
%
%     FIXED LEADER
%
%     n_ensembles - Number of ensemble is buffer
%     n_data - Number of data fields in buffer
%     firmver - Firmware version
%     firmrev - Firmware revision
%     sysconf - System configuration bytes
%     SymData - Flag for Real or Symulated data (0 for real data)
%     LagLength - Time period between sound pulses
%     usedbeams - Number of used beams
%     nbins - number of bins (1-128)
%     pingperens - pings per ensemble (0-16384)
%     binsize - bin size in cm (1-6400)
%     blnk - blanking in cm (0-9999)
%     watermode - Water profiling mode
%     minthrsh - Minimum threshold correlation in counts (1-256)
%     ncodrep - Code repetitions in transmit pulse in counts (1-256)
%     minpercgood - Minimum percentage of good pings in one ensemble (1-100)
%     maxerrvel - Maximum value of error velocity in mm/s (1-5000 mm/s)
%     Tbetweenpng -  %Time between two pings in cs (reading min, secs and cs)
%     corinfo - Coordinate information bytes
%     headalign - Head physical alignment correction in 0.01 degrees (-179.99 to 180.00)
%     headbias - Head magnetic bias correction in 0.01 degrees (-179.99 to 180.00)
%     sensource - Sensor source information (see manual for explanation)
%     senavail - Sensor availability info  (see manual for explanation)
%     distmidbin1 - Distance to middle of first bin in cm (0-65535)
%     lngthtranspulse - Length of the transmitted pulse in cm (0-65535)
%     watrefbins - Vector with begin and end bin for averaging to determine Water layer reference (1-128)
%     mintarget - Minimum for false target rejection in counts(0-255)
%     lowlattrig - Skip CX-command setting
%     distpulse - Distance between pulse repetitions in cm (0-65535) dependent on WM command
%     cpuserial - CPU serial number
%     bandwidth - Bandwidth (WB command)
%     syspower - System power in counts (CQ-command, only affects 75 and 150 KHz systems)
%     basefreqid - Base frequency index (only for Navigators)
%     serial - ADCP serial number (REMUS only)
%     beamangle - Beam angle, only for HADCP's
%
%     VARIABLE LEADER
%
%     ensnum - Ensemble number without rollover correction
%     BITcheck - Result of Built in test (see manual for explanation)
%     speedsound - Speed of sound in m/s (1400-1600)
%     depthtransd - Depth of transducer in decimeters (1-9999)
%     heading - Heading in 0.01 degrees (000.00-359.99)
%     pitch - Pitch in 0.01 degrees (-20.00 to 20.00)
%     roll - Roll in 0.01 degrees (-20.00 to 20.00)
%     salinity - Salinity in parts per thousands (0-40)
%     temperature - Temeperature in 0.01 degrees celcius (-5.00 to 40.00)
%     prepingT - Pre ping time in cs
%     headstd - Standard deviation in heading in degrees (0-180)
%     pitchstd - Standard deviation in pitch in 0.1 degrees (0.0-20.0)
%     rollstd - Standard deviation in roll in 0.1 degrees (0.0-20.0)
%     ADC - ADC channels, each columns is one channel
%     errorstat - Error status check (for explanation see manual)
%     pressure - Pressure in decapascal (0-4'294'967'295)
%     pressurevar - Variance in pressure in decapascal (0-4'294'967'295)
%     timeV - Gives date in Matlab time vector
%
%     ARRAY DATA
%
%     velocity
%     echo
%     correlation
%     percentage good
%
%     BOTTOM-TRACKING DATA
%
%     btpingperens - Bottom tracking pings per ensemble (0-999)
%     reacqdelay - Delay in number of ensembles before reacquiring (0-999)
%     mincormag                                %Minimum for correlation magnitude in counts (0-255)
%     minevampl                                %Minimum evaluation amplitude in counts (1-255)
%     btminpergood                             %Minimum percentage of good bt pings 
%     btmode                                   %BT mode
%     btmaxerrv                                %Maximum bt error velocity in mm/s (0-5000)
%     btrange                                  %bt range of beam 1,2,3,4 in cm (0-65535)
%     btvel                                    %bt velocty of beam 1,2,3,4 in mm/s (-32768 to 32768)
%     btcor                                    %bt correlation magnitude beam 1,2,3,4 in counts (0-255)
%     btevampl                                 %bt evaluation amplitude for strength of bottom echo beam 1,2,3,4 in counts (0-255)
%     btpercgood                               %bt percentage of good pings in beam 1,2,3,4 (0-100)
%     minlyrsize                               %Minimum size of ref layer in dm (0-9999)
%     nearbnd                                  %Near boundary of ref layer in dm (0-9999)
%     farbnd                                   %Far boundary of ref layer in dm (0-9999) 
%     reflyrvel                                %Reference layer velocity 1,2,3,4 in mm/s (-32768 to 32768)
%     reflyrcor                                %Reference layer correlation magnitude beam 1,2,3,4 in counts (0-255)
%     reflyrint                                %Reference layer echo intensity beam 1,2,3,4 in counts (0-255)
%     reflyrpergood                            %Reference layer percentage good pings in beam 1,2,3,4 (0-100)
%     maxdepth                                 %bt maximum depth in dm (80-9999)
%     rssiamp                                  %Received signal strength indicator in beam 1,2,3,4 in counts (0-255), 1 count ca. 0.45 dB
%     gain                                     %Gain level for shallow water
%
%     WINRIVER
%    
%     gga_siz
%     gga_dt
%     gga_header
%     gga_lat
%     gga_long
%     gga_qual
%     gga_nsats
%     gga_hdop
%     gga_altitude
%     gga_altunit
%     gga_geoid
%     gga_geoidunit
%     gga_agedgps
%     gga_refid
%
%
    
%
% PD0 Methods:
%
%

    
    
    properties
        buf
        n_ensembles
        max_nbins;
        max_nbeams;
    end

    properties(Dependent)
        
        %% fixed leader
        firmver;                                 %Firmware version
        firmrev;                                 %Firmware revision
        sysconf;                                 %System configuration (see manual for explanation)
                                                 % Does this display?       
        SymData;                                 %Flag for Real or Symulated data (0 for real data)
        LagLength;                               %Time period between sound pulses
        usedbeams;                               %Number of used beams
        nbins;                                   %number of bins (1-128)
        pingperens;                              %pings per ensemble (0-16384)
        binsize;                                 %bin size in cm (1-6400)
        blnk;                                    %blanking in cm (0-9999)
        watermode;                               %Water profiling mode
        minthrsh;                                %Minimum threshold correlation in counts (1-256)
        ncodrep;                                 %Code repetitions in transmit pulse in counts (1-256)
        minpercgood;                             %Minimum percentage of good pings in one ensemble (1-100)
        maxerrvel;                               %Maximum value of error velocity in mm/s (1-5000 mm/s)
        Tbetweenpng;                             %Time between two pings in cs (reading min, secs and cs)
        corinfo;                                 %Coordinate information (see manual for explanation)
        headalign;                               %Head physical alignment correction in 0.01 degrees (-179.99 to 180.00)
        headbias;                                %Head magnetic bias correction in 0.01 degrees (-179.99 to 180.00)
        sensource;                               %Sensor source information (see manual for explanation)
        senavail;                                %Sensor availability info  (see manual for explanation)
        distmidbin1;                             %Distance to middle of first bin in cm (0-65535)
        lngthtranspulse;                         %Length of the transmitted pulse in cm (0-65535)
        watrefbins;                              %Vector with begin and end bin for averaging to determine Water layer reference (1-128)
        mintarget;                               %Minimum for false target rejection in counts(0-255)
        lowlattrig;                              %Skip CX-command setting
        distpulse;                               %Distance between pulse repetitions in cm (0-65535) dependent on WM command
        cpuserial;                               %CPU serial number
        bandwidth;                               %Bandwidth (WB command)
        syspower;                                %System power in counts (CQ-command, only affects 75 and 150 KHz systems)
        basefreqid;                              %Base frequency index (only for Navigators)
        serial;                                  %ADCP serial number (REMUS only)
        beamangle;                               %Beam angle, only for HADCP's
        
        %% variable leader
        ensnum;                                  %Ensemble number without rollover correction
        BITcheck;                                %Result of Built in test (see manual for explanation)
        speedsound;                              %Speed of sound in m/s (1400-1600)
        depthtransd;                             %Depth of transducer in decimeters (1-9999)
        heading;                                 %Heading in 0.01 degrees (000.00-359.99)
        pitch;                                   %Pitch in 0.01 degrees (-20.00 to 20.00)
        roll;                                    %Roll in 0.01 degrees (-20.00 to 20.00)
        salinity;                                %Salinity in parts per thousands (0-40)
        temperature;                             %Temeperature in 0.01 degrees celcius (-5.00 to 40.00)
        prepingT;                                %Pre ping time in cs
        headstd;                                 %Standard deviation in heading in degrees (0-180)
        pitchstd;                                %Standard deviation in pitch in 0.1 degrees (0.0-20.0)
        rollstd;                                 %Standard deviation in roll in 0.1 degrees (0.0-20.0)
        ADC;                                     %ADC channels, each columns is one channel
        errorstat                                %Error status check (for explanation see manual)
        pressure;                                %Pressure in decapascal (0-4'294'967'295)
        pressurevar;                             %Variance in pressure in decapascal (0-4'294'967'295)
        timeV;                                   %Gives date in Matlab vector

        %% arrray data
        velocity;                                %Velocity data
        echo;                                    %Echo intensity counts
        percentage_good;                         %Percentage good data
        correlation;                             %Correlation data
        
        %% bottom tracking
        btpingperens                             %Bottom tracking pings per ensemble (0-999)
        reacqdelay                               %Delay in number of ensembles before reacquiring (0-999)
        mincormag                                %Minimum for correlation magnitude in counts (0-255)
        minevampl                                %Minimum evaluation amplitude in counts (1-255)
        btminpergood                             %Minimum percentage of good bt pings 
        btmode                                   %BT mode
        btmaxerrv                                %Maximum bt error velocity in mm/s (0-5000)
        btrange                                  %bt range of beam 1,2,3,4 in cm (0-65535)
        btvel                                    %bt velocty of beam 1,2,3,4 in mm/s (-32768 to 32768)
        btcor                                    %bt correlation magnitude beam 1,2,3,4 in counts (0-255)
        btevampl                                 %bt evaluation amplitude for strength of bottom echo beam 1,2,3,4 in counts (0-255)
        btpercgood                               %bt percentage of good pings in beam 1,2,3,4 (0-100)
        minlyrsize                               %Minimum size of ref layer in dm (0-9999)
        nearbnd                                  %Near boundary of ref layer in dm (0-9999)
        farbnd                                   %Far boundary of ref layer in dm (0-9999) 
        reflyrvel                                %Reference layer velocity 1,2,3,4 in mm/s (-32768 to 32768)
        reflyrcor                                %Reference layer correlation magnitude beam 1,2,3,4 in counts (0-255)
        reflyrint                                %Reference layer echo intensity beam 1,2,3,4 in counts (0-255)
        reflyrpergood                            %Reference layer percentage good pings in beam 1,2,3,4 (0-100)
        maxdepth                                 %bt maximum depth in dm (80-9999)
        rssiamp                                  %Received signal strength indicator in beam 1,2,3,4 in counts (0-255), 1 count ca. 0.45 dB
        gain                                     %Gain level for shallow water

        
        %% winriver data

        gga_siz
        gga_dt
        gga_header
        gga_lat
        gga_long
        gga_qual
        gga_nsats
        gga_hdop
        gga_altitude
        gga_altunit
        gga_geoid
        gga_geoidunit
        gga_agedgps
        gga_refid
        
    end
    properties(Access=private)
        ens_pos
        ens_bytes
        ens_ndat
        data_ensid
        data_offset
        data_headers
        array_subs
        array_idx
        n_data
    end

    
    
    methods
        %% CONSTRUCTOR
        function obj=PD0(varargin)
            if nargin>0
                if isa(varargin{1},'uint8') && iscolumn(varargin{1}) % Construct with a buffer
                    obj.buf=varargin{1};
                elseif ischar(varargin{1}) && isrow(varargin{1}) && exist(varargin{1},'file') % Construct with a filename
                    [fid, message]=fopen(varargin{1},'r','l');
                    if fid==-1
                        error(message);
                    end
                    obj.buf=fread(fid,'*uint8');
                    fclose(fid);
                elseif isa(varargin{1},'rdi.PD0') % Copy constructor (copies handle)
                    obj=varargin{1};
                else
                    warning('rdi:PD0:WrongConstruction','Did not understand construction argument, constructing empty')
                end
            end
        end
        
        %% SETTERS AND GETTERS
        function set.buf(obj,val)
            assert(isa(val,'uint8'));
            assert(iscolumn(val));
            obj.buf=val;
            obj.init();
        end
        
        %% Fixed leader
        function val=get.firmver(obj),          val=obj.get_scalar_data(rdi.headers.fixed_leader,2,'uint8'); end
        function val=get.firmrev(obj),          val=obj.get_scalar_data(rdi.headers.fixed_leader,3,'uint8'); end
        function val=get.sysconf(obj),          val=dec2bin(obj.get_scalar_data(rdi.headers.fixed_leader,4,'uint16'),16)'=='1'; end
        function val=get.SymData(obj),          val=obj.get_scalar_data(rdi.headers.fixed_leader,6,'uint8');end
        function val=get.LagLength(obj),        val=obj.get_scalar_data(rdi.headers.fixed_leader,7,'uint8');end
        function val=get.usedbeams(obj),        val=obj.get_scalar_data(rdi.headers.fixed_leader,8,'uint8');end
        function val=get.nbins(obj),            val=obj.get_scalar_data(rdi.headers.fixed_leader,9,'uint8');end
        function val=get.pingperens(obj),       val=obj.get_scalar_data(rdi.headers.fixed_leader,10,'uint16'); end
        function val=get.binsize(obj),          val=obj.get_scalar_data(rdi.headers.fixed_leader,12,'uint16'); end
        function val=get.blnk(obj),             val=obj.get_scalar_data(rdi.headers.fixed_leader,14,'uint16'); end
        function val=get.watermode(obj),        val=obj.get_scalar_data(rdi.headers.fixed_leader,16,'uint8'); end        
        function val=get.minthrsh(obj),         val=obj.get_scalar_data(rdi.headers.fixed_leader,17,'uint8'); end
        function val=get.ncodrep(obj),          val=obj.get_scalar_data(rdi.headers.fixed_leader,18,'uint8'); end
        function val=get.minpercgood(obj),      val=obj.get_scalar_data(rdi.headers.fixed_leader,19,'uint8'); end
        function val=get.maxerrvel(obj),        val=obj.get_scalar_data(rdi.headers.fixed_leader,20,'uint16'); end
        function val=get.Tbetweenpng(obj),      val=uint16(obj.get_scalar_data(rdi.headers.fixed_leader,22,'uint8'))*6000+...        
                                                    uint16(obj.get_scalar_data(rdi.headers.fixed_leader,23,'uint8'))*100+...        
                                                    uint16(obj.get_scalar_data(rdi.headers.fixed_leader,24,'uint8')); end
        function val=get.corinfo(obj),          val=dec2bin(obj.get_scalar_data(rdi.headers.fixed_leader,25,'uint8'),8)'=='1'; end
        function val=get.headalign(obj),        val=obj.get_scalar_data(rdi.headers.fixed_leader,26,'int16'); end
        function val=get.headbias(obj),         val=obj.get_scalar_data(rdi.headers.fixed_leader,28,'int16'); end
        function val=get.sensource(obj),        val=dec2bin(obj.get_scalar_data(rdi.headers.fixed_leader,30,'uint8'),8)'=='1'; end
        function val=get.senavail(obj),         val=obj.get_scalar_data(rdi.headers.fixed_leader,31,'uint8'); end
        function val=get.distmidbin1(obj),      val=obj.get_scalar_data(rdi.headers.fixed_leader,32,'uint16'); end
        function val=get.lngthtranspulse(obj),  val=obj.get_scalar_data(rdi.headers.fixed_leader,34,'uint16'); end
        function val=get.watrefbins(obj),       val=obj.get_scalar_data(rdi.headers.fixed_leader,36,'uint16'); end
        function val=get.mintarget(obj),        val=obj.get_scalar_data(rdi.headers.fixed_leader,38,'uint8'); end
        function val=get.lowlattrig(obj),       val=obj.get_scalar_data(rdi.headers.fixed_leader,39,'uint8'); end
        function val=get.distpulse(obj),        val=obj.get_scalar_data(rdi.headers.fixed_leader,40,'uint16'); end
        function val=get.cpuserial(obj),        val=[obj.get_scalar_data(rdi.headers.fixed_leader,42,'uint8');...       
                                                    obj.get_scalar_data(rdi.headers.fixed_leader,43,'uint8');...       
                                                    obj.get_scalar_data(rdi.headers.fixed_leader,44,'uint8');...        
                                                    obj.get_scalar_data(rdi.headers.fixed_leader,45,'uint8');...        
                                                    obj.get_scalar_data(rdi.headers.fixed_leader,46,'uint8');...        
                                                    obj.get_scalar_data(rdi.headers.fixed_leader,47,'uint8');...        
                                                    obj.get_scalar_data(rdi.headers.fixed_leader,48,'uint8');...        
                                                    obj.get_scalar_data(rdi.headers.fixed_leader,49,'uint8')]; end
        function val=get.bandwidth(obj),        val=obj.get_scalar_data(rdi.headers.fixed_leader,50,'uint16'); end
        function val=get.syspower(obj),         val=obj.get_scalar_data(rdi.headers.fixed_leader,52,'uint8'); end
        function val=get.basefreqid(obj),       val=obj.get_scalar_data(rdi.headers.fixed_leader,53,'uint8'); end
        function val=get.serial(obj),           val=[obj.get_scalar_data(rdi.headers.fixed_leader,54,'uint8');...
                                                     obj.get_scalar_data(rdi.headers.fixed_leader,55,'uint8');...
                                                     obj.get_scalar_data(rdi.headers.fixed_leader,56,'uint8');...
                                                     obj.get_scalar_data(rdi.headers.fixed_leader,57,'uint8')]; end
        function val=get.beamangle(obj),        val=obj.get_scalar_data(rdi.headers.fixed_leader,58,'uint8'); end
        
        %% Variable leader
        function val=get.ensnum(obj),           val=uint32(obj.get_scalar_data(rdi.headers.variable_leader,2,'uint16'))+...        
                                                    uint32(obj.get_scalar_data(rdi.headers.variable_leader,11,'uint8'))*65535; end
        function val=get.BITcheck(obj),         val=dec2bin(obj.get_scalar_data(rdi.headers.variable_leader,12,'uint16'),16)'=='1'; end
        function val=get.speedsound(obj),       val=obj.get_scalar_data(rdi.headers.variable_leader,14,'uint16'); end
        function val=get.depthtransd(obj),      val=obj.get_scalar_data(rdi.headers.variable_leader,16,'uint16'); end
        function val=get.heading(obj),          val=obj.get_scalar_data(rdi.headers.variable_leader,18,'uint16'); end
        function val=get.pitch(obj),            val=obj.get_scalar_data(rdi.headers.variable_leader,20,'int16'); end
        function val=get.roll(obj),             val=obj.get_scalar_data(rdi.headers.variable_leader,22,'int16'); end
        function val=get.salinity(obj),         val=obj.get_scalar_data(rdi.headers.variable_leader,24,'uint16'); end
        function val=get.temperature(obj),      val=obj.get_scalar_data(rdi.headers.variable_leader,26,'int16'); end
        function val=get.prepingT(obj),         val=uint16(obj.get_scalar_data(rdi.headers.variable_leader,28,'uint8'))*600+...
                                                    uint16(obj.get_scalar_data(rdi.headers.variable_leader,29,'uint8'))*100+...
                                                    uint16(obj.get_scalar_data(rdi.headers.variable_leader,30,'uint8')); end
        function val=get.headstd(obj),          val=obj.get_scalar_data(rdi.headers.variable_leader,31,'uint8'); end
        function val=get.pitchstd(obj),         val=obj.get_scalar_data(rdi.headers.variable_leader,32,'uint8'); end
        function val=get.rollstd(obj),          val=obj.get_scalar_data(rdi.headers.variable_leader,33,'uint8'); end
        function val=get.ADC(obj),              val=[obj.get_scalar_data(rdi.headers.variable_leader,34,'uint8');...
                                                     obj.get_scalar_data(rdi.headers.variable_leader,35,'uint8');...
                                                     obj.get_scalar_data(rdi.headers.variable_leader,36,'uint8');...
                                                     obj.get_scalar_data(rdi.headers.variable_leader,37,'uint8');...
                                                     obj.get_scalar_data(rdi.headers.variable_leader,38,'uint8');...
                                                     obj.get_scalar_data(rdi.headers.variable_leader,39,'uint8');...
                                                     obj.get_scalar_data(rdi.headers.variable_leader,40,'uint8');...
                                                     obj.get_scalar_data(rdi.headers.variable_leader,41,'uint8')]; end
        function val=get.errorstat(obj),        val=[dec2bin(obj.get_scalar_data(rdi.headers.variable_leader,42,'uint8'),8)'=='1';...
                                                     dec2bin(obj.get_scalar_data(rdi.headers.variable_leader,43,'uint8'),8)'=='1';...
                                                     dec2bin(obj.get_scalar_data(rdi.headers.variable_leader,44,'uint8'),8)'=='1';...
                                                     dec2bin(obj.get_scalar_data(rdi.headers.variable_leader,45,'uint8'),8)'=='1']; end
        function val=get.pressure(obj),         val=obj.get_scalar_data(rdi.headers.variable_leader,48,'uint32'); end
        function val=get.pressurevar(obj),      val=obj.get_scalar_data(rdi.headers.variable_leader,52,'uint32'); end
        function val=get.timeV(obj),            val=[double(obj.get_scalar_data(rdi.headers.variable_leader,57,'uint8'))*100+...        
                                                     double(obj.get_scalar_data(rdi.headers.variable_leader,58,'uint8'));...
                                                     double(obj.get_scalar_data(rdi.headers.variable_leader,59,'uint8'));...
                                                     double(obj.get_scalar_data(rdi.headers.variable_leader,60,'uint8'));...
                                                     double(obj.get_scalar_data(rdi.headers.variable_leader,61,'uint8'));...
                                                     double(obj.get_scalar_data(rdi.headers.variable_leader,62,'uint8'));...
                                                     double(obj.get_scalar_data(rdi.headers.variable_leader,63,'uint8'))+...
                                                     double(obj.get_scalar_data(rdi.headers.variable_leader,64,'uint8'))/100]; end

        %% array data
        function val=get.velocity(obj),         val=obj.get_array_data('int16',rdi.headers.velocity); end
        function val=get.echo(obj),             val=obj.get_array_data('uint8',rdi.headers.echo); end
        function val=get.percentage_good(obj),  val=obj.get_array_data('uint8',rdi.headers.percentage_good); end
        function val=get.correlation(obj),      val=obj.get_array_data('uint8',rdi.headers.correlation); end
        
        %% Bottom tracking data 
        function val=get.btpingperens(obj),     val=obj.get_scalar_data(rdi.headers.bottom_tracking,2,'uint16'); end
        function val=get.reacqdelay(obj),       val=obj.get_scalar_data(rdi.headers.bottom_tracking,4,'uint16'); end
        function val=get.mincormag(obj),        val=obj.get_scalar_data(rdi.headers.bottom_tracking,6,'uint8'); end
        function val=get.minevampl(obj),        val=obj.get_scalar_data(rdi.headers.bottom_tracking,7,'uint8'); end
        function val=get.btminpergood(obj),     val=obj.get_scalar_data(rdi.headers.bottom_tracking,8,'uint8'); end
        function val=get.btmode(obj),           val=obj.get_scalar_data(rdi.headers.bottom_tracking,9,'uint8'); end
        function val=get.btmaxerrv(obj),        val=obj.get_scalar_data(rdi.headers.bottom_tracking,10,'uint16'); end
        function val=get.btrange(obj),          val=uint32(cat(3,obj.get_scalar_data(rdi.headers.bottom_tracking,16,'uint16'),...
                                                                 obj.get_scalar_data(rdi.headers.bottom_tracking,18,'uint16'),...
                                                                 obj.get_scalar_data(rdi.headers.bottom_tracking,20,'uint16'),...
                                                                 obj.get_scalar_data(rdi.headers.bottom_tracking,22,'uint16')))+...
                                                    uint32(cat(3,obj.get_scalar_data(rdi.headers.bottom_tracking,77,'uint8'),...
                                                                 obj.get_scalar_data(rdi.headers.bottom_tracking,78,'uint8'),...
                                                                 obj.get_scalar_data(rdi.headers.bottom_tracking,79,'uint8'),...
                                                                 obj.get_scalar_data(rdi.headers.bottom_tracking,80,'uint8')))*65536;
            if isempty(val), val=[];end
        end
        function val=get.btvel(obj),            val=cat(3,obj.get_scalar_data(rdi.headers.bottom_tracking,24,'int16'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,26,'int16'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,28,'int16'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,30,'int16'));
            if isempty(val), val=[];end
        end
        function val=get.btcor(obj),            val=cat(3,obj.get_scalar_data(rdi.headers.bottom_tracking,32,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,33,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,34,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,35,'uint8'));
            if isempty(val), val=[];end
        end
        function val=get.btevampl(obj),         val=cat(3,obj.get_scalar_data(rdi.headers.bottom_tracking,36,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,37,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,38,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,39,'uint8'));
            if isempty(val), val=[];end
        end
        function val=get.btpercgood(obj),       val=cat(3,obj.get_scalar_data(rdi.headers.bottom_tracking,40,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,41,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,42,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,43,'uint8'));
            if isempty(val), val=[];end
        end
        function val=get.minlyrsize(obj),       val=obj.get_scalar_data(rdi.headers.bottom_tracking,44,'uint16'); end
        function val=get.nearbnd(obj),          val=obj.get_scalar_data(rdi.headers.bottom_tracking,46,'uint16'); end
        function val=get.farbnd(obj),           val=obj.get_scalar_data(rdi.headers.bottom_tracking,48,'uint16'); end
        function val=get.reflyrvel(obj),        val=cat(3,obj.get_scalar_data(rdi.headers.bottom_tracking,50,'int16'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,52,'int16'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,54,'int16'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,56,'int16'));
            if isempty(val), val=[];end
        end
        function val=get.reflyrcor(obj),        val=cat(3,obj.get_scalar_data(rdi.headers.bottom_tracking,58,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,59,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,60,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,61,'uint8'));
            if isempty(val), val=[];end
        end
        function val=get.reflyrint(obj),        val=cat(3,obj.get_scalar_data(rdi.headers.bottom_tracking,62,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,63,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,64,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,65,'uint8'));
            if isempty(val), val=[];end
        end
        function val=get.reflyrpergood(obj),    val=cat(3,obj.get_scalar_data(rdi.headers.bottom_tracking,66,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,67,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,68,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,69,'uint8'));
            if isempty(val), val=[];end
        end
        function val=get.maxdepth(obj),         val=obj.get_scalar_data(rdi.headers.bottom_tracking,70,'uint16'); end
        function val=get.rssiamp(obj),          val=cat(3,obj.get_scalar_data(rdi.headers.bottom_tracking,72,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,73,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,74,'uint8'),...
                                                          obj.get_scalar_data(rdi.headers.bottom_tracking,75,'uint8'));
            if isempty(val), val=[];end
        end
        function val=get.gain(obj),             val=obj.get_scalar_data(rdi.headers.bottom_tracking,76,'uint8'); end

        %% winriver data (we assume either wr1 or wr2_v1 or wr2_v2 data. Nothing simultaneous in one ensemble
        
        function val=get.gga_siz(obj)
            val=obj.get_scalar_data(rdi.headers.NMEA_Specific_GGA_v2,5,'uint16');
        end
        function val=get.gga_dt(obj)
            val=obj.get_scalar_data(rdi.headers.NMEA_Specific_GGA_v2,6,'double');
        end
        function val=get.gga_header(obj)
            val=char(obj.get_scalar_data(rdi.headers.NMEA_Specific_GGA_v2,14,'uint8'));
        end
%         function val=get.gga_lat
%         function val=get.gga_long
%         function val=get.gga_qual
%         function val=get.gga_nsats
%         function val=get.gga_hdop
%         function val=get.gga_altitude
%         function val=get.gga_altunit
%         function val=get.gga_geoid
%         function val=get.gga_geoidunit
%         function val=get.gga_agedgps
%         function val=get.gga_refid
    end
    
    %% PRIVATE METHODS
    methods(Access=private)
        init(obj);  
        out=parse_blocks(obj,pos,type)
        function data_idx=find_data(obj,type)
            data_idx=nan(1,obj.n_ensembles);
            ftype=find(obj.data_headers==type);
            data_idx(obj.data_ensid(ftype))=ftype; % Keeps last, if multiple data blocks of same type. nan means data not found in ensemble
        end

        function val=get_array_data(obj,type,header) 
            nullval=obj.nullval(type);
            val=nullval(ones(obj.max_nbins,obj.n_ensembles,obj.max_nbeams)); % Initialize
            if isempty(val), return, end;
            
            ndat=double(obj.nbins).*double(obj.usedbeams);
            ensidx=obj.array_subs(2,:);
            idces=obj.array_idx;
            
            % Velocity position
            pos=obj.find_data(header); % Locate all velocity data
                % account for possible ensembles without data (does this ever happen?)
                f_bad_ens=find(~isfinite(pos));
                [~,f_bad_idx,~]=intersect(ensidx,f_bad_ens); % Get index of all velocities which are missing
                ensidx(f_bad_idx)=[];
                ndat(f_bad_ens)=[];
            pos=obj.data_offset(pos(ensidx))+cumsum(obj.sizeof(type)*ones(1,sum(ndat)));
            cnvels=cumsum([0 ndat]);
            pos=pos-cnvels(ensidx)*obj.sizeof(type);
            idces(f_bad_idx)=[];
            val(idces)=obj.parse_blocks(pos,type);
        end
        
        function val=get_scalar_data(obj,header,idx,type)
            dataidx=obj.find_data(header);
            fgood=isfinite(dataidx);
            if ~any(fgood), val=[]; return, end;
            nullval=obj.nullval(type);
            val=nullval(ones(1,obj.n_ensembles));
            val(obj.data_ensid(dataidx(fgood)))=obj.parse_blocks(obj.data_offset(dataidx(fgood))+idx,type);  
        end
    end
    methods(Access=private, Static)
        function nb=sizeof(type)
            nb=cast(1,type); %#ok<NASGU>
            nb=whos('nb');
            nb=nb.bytes;
        end
        function val=nullval(type)
            if strcmp(type(1:3),'int'), val=intmin(type); 
            elseif strcmp(type(1:4),'uint'), val = intmax(type);
            elseif strcmp(type,{'single','double'}), val=nan;
            else val=0;
            end
        end
    end
end