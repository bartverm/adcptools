classdef ADCP_Type
% Enumeration defining different ADCP models
    enumeration
        Unknown
        ChannelMaster
        ExplorerPhasedArray
        ExplorerPiston
        SentinelV
        MonitorV
        RioGrande
        RiverRay
        StreamPro
        
        % Workhorse below:
        Sentinel
        Mariner
        Monitor
        QuarterMaster1500
        QuarterMaster3000
        QuarterMaster6000
        QuarterMaster1500ModBeams
        LongRanger75
        LongRanger1500
        LongRanger3000
    end
end