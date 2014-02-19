classdef headers < uint16
    enumeration
        velocity (256)
        echo (768)
        correlation (512)
        percentage_good (1024)
        bottom_tracking (1536)
        variable_leader (128)
        fixed_leader (0)
        WinRiverII (8226)
        NMEA_DBT (8448)
        NMEA_GGA (8449)
        NMEA_VTG (8450)
        NMEA_GSA (8451)
        NMEA_HDT (8452)
        Bad_header (65535)
    end
end