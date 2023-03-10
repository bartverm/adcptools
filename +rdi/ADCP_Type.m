classdef ADCP_Type < uint8
% Enumeration defining different ADCP models
%

% This information comes from the PDDecoder software from Teledyne Marine
% which used to be on github for a while
    enumeration
        Unknown(0),
        BROADBAND_5(5),                            % Old broadband systems
        WORKHORSE_8(8),                            % Old Workhorse
        NAVIGATOR_9(9),                            % Navigator
        RIO_GRANDE_10(10),                         % RioGrande
        WH_HORIZONTAL_11(11),                      % Workhorse horizontal system
        OCEAN_SURVEYOR_14(14),                     % Ocean Surveyor
        WORKHORSE_16(16),                          % Workhorse
        NAVIGATOR_19(19),                          % Navigator
        OCEAN_SURVEYOR_23(23),                      % Ocean Surveyor
        CHANNELMASTER_28(28),                      % ChannelMaster
        STREAMPRO_31(31),                          % StreamPro
        EXPLORER_34(34),                           % Explorer
        NAVIGATOR_37(37),                          % Navigator
        DVS_41(41),                                % DVS
        WORKHORSE_43(43),                          % Workhorse
        RIVERRAY_44(44),                           % RiverRay
        SENTINELV_47(47),                          % SentinelV
        WORKHORSE_50(50),                          % Workhorse
        WORKHORSE_51(51),                          % Workhorse
        WORKHORSE_52(52),                          % Workhorse
        NAVIGATOR_53(53),                          % Workhorse
        DVS_55(55),                                % DVS
        RIVERPRO_56(56),                           % RiverPro(5 beams) or RioPro(4 beams)
        MERIDIAN_59(59),                           % Meridian
        PINNACLE_61(61),                           % Pinnacle
        SENTINELV_66(66),                          % SentinelV RT
        PATHFINDER_67(67),                         % Pathfinder
        PIONEER_73(73),                            % Pioneer
        TASMAN_74(74),                             % TASMAN
        WAYFINDER_76(76),                          % Wayfinder
        WORKHORSE_77(77),                          % Workhorse II
        WORKHORSE_78(78),                          % Workhorse II
    end
end