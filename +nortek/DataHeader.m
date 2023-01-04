classdef DataHeader < uint8
    enumeration
        Unknown (0x00)
        ConfigurationString(0x10)
        AD2CP (0xA5)
        Burst (0x15)
        Average (0x16)
        BottomTracking (0x17)
        Interleaved (0x18)
        BurstAltimeter (0x1A)
        DVLBottomTracking (0x1B)
        EchoSounder (0x1C)
        DVLWaterTracking (0x1D)
        Altimeter (0x1E)
        AverageAltimeter (0x1F)
        String (0xA0)
        
        % PTP base64 strings
        PTP0 (0xB0)
        PTP1 (0xB1)
        PTP2 (0xB2)
        PTP3 (0xB3)
        PTP4 (0xB4)
        PTP5 (0xB5)
        PTP6 (0xB6)
        PTP7 (0xB7)
        PTP8 (0xB8)
        PTP9 (0xB9)
        PTPa (0xBA)
        PTPb (0xBB)
        PTPc (0xBC)
        PTPd (0xBD)
        PTPe (0xBE)
        PTPf (0xBF)
    end
end