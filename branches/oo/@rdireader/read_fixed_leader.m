function read_fixed_leader(obj,fpos)
fseek(obj.fid,fpos,-1);
LeadID=fread(obj.fid,2,'*uint8');
if any(LeadID~=[0;0])
    error('rdireader:read_fixed_leader:wrongID','Fixed Leader ID seems to be wrong')
end
obj.adcp.manufacturer.firmver=fread(obj.fid,1,'*uint8');                                 %Firmware version
obj.adcp.manufacturer.firmrev=fread(obj.fid,1,'*uint8');                                 %Firmware revision
obj.adcp.manufacturer.sysconf=num2str(fread(obj.fid,16,'*ubit1'))';                      %System configuration (see manual for explanation)
% obj.adcp.manufacturer.sysconfstr=sysinterp(obj.adcp.manufacturer.sysconf);                     %Interprete system configuration
obj.adcp.manufacturer.SymData=fread(obj.fid,1,'*uint8');                                 %Flag for Real or Symulated data (0 for real data)
obj.adcp.manufacturer.LagLength=fread(obj.fid,1,'*uint8');                               %Time period between sound pulses
obj.adcp.manufacturer.usedbeams=fread(obj.fid,1,'*uint8');                               %Number of used beams
obj.adcp.manufacturer.nbins=fread(obj.fid,1,'*uint8');                                   %number of bins (1-128)
obj.adcp.manufacturer.pingperens=fread(obj.fid,1,'*uint16');                             %pings per ensemble (0-16384)
obj.adcp.manufacturer.binsize=fread(obj.fid,1,'*uint16');                                %bin size in cm (1-6400)
obj.adcp.manufacturer.blnk=fread(obj.fid,1,'*uint16');                                   %blanking in cm (0-9999)
fseek(obj.fid,1,0);                                                            %skip data processing mode (always one)
obj.adcp.manufacturer.minthrsh=fread(obj.fid,1,'*uint8');                                %Minimum threshold correlation in counts (1-256)
obj.adcp.manufacturer.ncodrep=fread(obj.fid,1,'*uint8');                                 %Code repetitions in transmit pulse in counts (1-256)
obj.adcp.manufacturer.minpercgood=fread(obj.fid,1,'*uint8');                             %Minimum percentage of good pings in one ensemble (1-100)
obj.adcp.manufacturer.maxerrvel=fread(obj.fid,1,'*uint16');                              %Maximum value of error velocity in mm/s (1-5000 mm/s)
obj.adcp.manufacturer.Tbetweenpng=fread(obj.fid,1,'uint8=>uint16')*6000+...
    fread(obj.fid,1,'uint8=>uint16')*100+fread(obj.fid,1,'uint8=>uint16');         %Time between two pings in cs (reading min, secs and cs)
obj.adcp.manufacturer.corinfo=num2str(fread(obj.fid,8,'*ubit1'))';                       %Coordinate information (see manual for explanation)
% obj.adcp.manufacturer.corstr{filenum}=corinterpr(obj.adcp.manufacturer.corinfo);                        %Interprete coordinate information
obj.adcp.manufacturer.headalign=fread(obj.fid,1,'*int16');                               %Head physical alignment correction in 0.01 degrees (-179.99 to 180.00)
obj.adcp.manufacturer.headbias=fread(obj.fid,1,'*int16');                                %Head magnetic bias correction in 0.01 degrees (-179.99 to 180.00)
obj.adcp.manufacturer.sensource=num2str(fread(obj.fid,8,'*ubit1'))';                     %Sensor source information (see manual for explanation)
obj.adcp.manufacturer.senavail=num2str(fread(obj.fid,8,'*ubit1'))';                      %Sensor availability info  (see manual for explanation)
obj.adcp.manufacturer.distmidbin1=fread(obj.fid,1,'*uint16');                            %Distance to middle of first bin in cm (0-65535)
obj.adcp.manufacturer.lngthtranspulse=fread(obj.fid,1,'*uint16');                        %Length of the transmitted pulse in cm (0-65535)
obj.adcp.manufacturer.watrefbins=fread(obj.fid,2,'*uint8')';                              %Vector with begin and end bin for averaging to determine Water layer reference (1-128)
obj.adcp.manufacturer.mintarget=fread(obj.fid,1,'*uint8');                               %Minimum for false target rejection in counts(0-255)
obj.adcp.manufacturer.lowlattrig=fread(obj.fid,1,'*uint8');                              %Skip CX-command setting
obj.adcp.manufacturer.distpulse=fread(obj.fid,1,'*uint16');                              %Distance between pulse repetitions in cm (0-65535) dependent on WM command
obj.adcp.manufacturer.cpuserial=fread(obj.fid,8,'*uint8')';                               %CPU serial number
obj.adcp.manufacturer.bandwidth=fread(obj.fid,1,'*uint16');                              %Bandwidth (WB command)
obj.adcp.manufacturer.syspower=fread(obj.fid,1,'*uint8');                                %System power in counts (CQ-command, only affects 75 and 150 KHz systems)
obj.adcp.manufacturer.basefreqid=fread(obj.fid,1,'*uint8');                              %Base frequency index (only for Navigators)
obj.adcp.manufacturer.serial=fread(obj.fid,4,'*uint8')';                                  %ADCP serial number (REMUS only)
obj.adcp.manufacturer.HADCPbeamangle=fread(obj.fid,1,'*uint8');                          %Beam angle, only for HADCP's
