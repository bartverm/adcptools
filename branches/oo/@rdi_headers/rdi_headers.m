classdef rdi_headers < uint16
    enumeration
        Velocity (256)
        Echo (768)
        Correlation (512)
        Percentage_good (1024)
        Bottom_tracking (1536)
        Variable_leader (128)
        WinRiverII (8226)
        NMEA_DBT (8448)
        NMEA_GGA (8449)
        NMEA_VTG (8450)
        NMEA_HDT (8452)
        Bad_header (65535)
    end
end