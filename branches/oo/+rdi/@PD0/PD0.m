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

%     PD0 Methods:

    
    
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
        
        %% external data
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
        function val=get.firmver(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+2);        
        end
        function val=get.firmrev(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+3);        
        end
        function val=get.sysconf(obj)
            val=dec2bin(obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+4,'uint16'),16);        
        end
        function val=get.SymData(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+6);        
        end
        function val=get.LagLength(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+7);        
        end
        function val=get.usedbeams(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+8);        
        end
        function val=get.nbins(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+9);        
        end
        function val=get.pingperens(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+10,'uint16');        
        end
        function val=get.binsize(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+12,'uint16');        
        end
        function val=get.blnk(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+14,'uint16');        
        end
        function val=get.watermode(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+16);        
        end        
        function val=get.minthrsh(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+17);        
        end
        function val=get.ncodrep(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+18);        
        end
        function val=get.minpercgood(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+19);        
        end
        function val=get.maxerrvel(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+20,'uint16');        
        end
        function val=get.Tbetweenpng(obj)
            val=uint16(obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+22))*6000+...        
                uint16(obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+23))*100+...        
                uint16(obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+24));        
        end
        function val=get.corinfo(obj)
            val=dec2bin(obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+25),8);        
        end
        function val=get.headalign(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+26,'int16');        
        end
        function val=get.headbias(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+28,'int16');        
        end
        function val=get.sensource(obj)
            val=dec2bin(obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+30),8);        
        end
        function val=get.senavail(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+31);        
        end
        function val=get.distmidbin1(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+32,'uint16');        
        end
        function val=get.lngthtranspulse(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+34,'uint16');        
        end
        function val=get.watrefbins(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+36,'uint16');        
        end
        function val=get.mintarget(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+38);        
        end
        function val=get.lowlattrig(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+39);        
        end
        function val=get.distpulse(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+40,'uint16');        
        end
        function val=get.cpuserial(obj)
            val=[obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+42),...       
                 obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+43),...       
            	 obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+44),...        
            	 obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+45),...        
           	     obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+46),...        
            	 obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+47),...        
            	 obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+48),...        
            	 obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+49)];      
        end
        function val=get.bandwidth(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+50,'uint16');        
        end
        function val=get.syspower(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+52);        
        end
        function val=get.basefreqid(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+53);        
        end
        function val=get.serial(obj)
            val=[obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+54),...
            	 obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+55),...
            	 obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+56),...
            	 obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+57)];
        end
        function val=get.beamangle(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.fixed_leader))+58);        
        end
        
        %% Variable leader
        function val=get.ensnum(obj)
            val=uint32(obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+2,'uint16'))+...        
                uint32(obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+11))*65535;        
        end
        function val=get.BITcheck(obj)
            val=dec2bin(obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+12,'uint16'),16);        
        end
        function val=get.speedsound(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+14,'uint16');        
        end
        function val=get.depthtransd(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+16,'uint16');        
        end
        function val=get.heading(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+18,'uint16');        
        end
        function val=get.pitch(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+20,'int16');        
        end
        function val=get.roll(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+22,'int16');        
        end
        function val=get.salinity(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+24,'uint16');        
        end
        function val=get.temperature(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+26,'int16');        
        end
        function val=get.prepingT(obj)
            val=uint16(obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+28))*600+...
                uint16(obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+29))*100+...
                uint16(obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+30));        
        end
        function val=get.headstd(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+31);        
        end
        function val=get.pitchstd(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+32);        
        end
        function val=get.rollstd(obj)
            val=obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+33);        
        end
        function val=get.ADC(obj)
            val=[obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+34),...
            	 obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+35),...
            	 obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+36),...
            	 obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+37),...
            	 obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+38),...
            	 obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+39),...
            	 obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+40),...
            	 obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+41)];
        end
        function val=get.errorstat(obj)
            val=[dec2bin(obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+42),8),...
            	 dec2bin(obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+43),8),...
            	 dec2bin(obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+44),8),...
            	 dec2bin(obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+45),8)];
        end
        function val=get.pressure(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+48,'uint32');        
        end
        function val=get.pressurevar(obj)
            val=obj.parse_blocks(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+52,'uint32');        
        end
        function val=get.timeV(obj)
            val=[double(obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+57))*100+...        
                 double(obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+58)),...
                 double(obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+59)),...
                 double(obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+60)),...
                 double(obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+61)),...
                 double(obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+62)),...
                 double(obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+63))+...
                 double(obj.buf(obj.data_offset(obj.find_data(rdi.headers.variable_leader))+64))/100];        
        end

        %% array data
        function val=get.velocity(obj)
            val=obj.get_array_data('int16',rdi.headers.velocity);
        end
        function val=get.echo(obj)
            val=obj.get_array_data('uint8',rdi.headers.echo);
        end
        function val=get.percentage_good(obj)
            val=obj.get_array_data('uint8',rdi.headers.percentage_good);
        end
        function val=get.correlation(obj)
            val=obj.get_array_data('uint8',rdi.headers.correlation);
        end
        
        
    end
    
    %% PRIVATE METHODS
    methods(Access=private)
        init(obj);  
        out=parse_blocks(obj,pos,type)
        function data_idx=find_data(obj,type)
            data_idx=nan(obj.n_ensembles,1);
            ftype=find(obj.data_headers==type);
            data_idx(obj.data_ensid(ftype))=ftype; % Keeps last, if multiple data blocks of same type. nan means data not found in ensemble
        end

        function val=get_array_data(obj,type,header) 
            if type(1)=='i', nullval=intmin(type); else nullval = intmax(type); end
            
            val=nullval(ones(obj.max_nbins,obj.n_ensembles,obj.max_nbeams)); % Initialize

            ndat=double(obj.nbins).*double(obj.usedbeams);
            ensidx=obj.array_subs(:,2);
            idces=obj.array_idx;
            
            % Velocity position
            pos=obj.find_data(header); % Locate all velocity data
                % account for possible ensembles without data (does this ever happen?)
                f_bad_ens=find(~isfinite(pos));
                [~,f_bad_idx,~]=intersect(ensidx,f_bad_ens); % Get index of all velocities which are missing
                ensidx(f_bad_idx)=[];
                ndat(f_bad_ens)=[];
            pos=obj.data_offset(pos(ensidx))+cumsum(rdi.PD0.sizeof(type)*ones(sum(ndat),1));
            cnvels=cumsum([0;ndat]);
            pos=pos-cnvels(ensidx)*rdi.PD0.sizeof(type);
            idces(f_bad_idx)=[];
            val(idces)=obj.parse_blocks(pos,type);
        end
    end
    methods(Access=private, Static)
        function nb=sizeof(type)
            nb=cast(1,type); %#ok<NASGU>
            nb=whos('nb');
            nb=nb.bytes;
        end
    end
end