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
        Unknown2 (0xB0) % Need to figure this out!
    end
end