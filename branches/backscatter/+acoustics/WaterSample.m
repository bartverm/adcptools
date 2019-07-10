classdef WaterSample < handle
    properties
        distance(1,1) double {mustBeFinite, mustBeNonnegative} = 0
        concentration(1,1) double {mustBeFinite, mustBeNonnegative} = 0
        time(1,1) datetime = datetime()
        distribution(1,1) acoustics.GrainSizeDistribution = acoustics.GrainSizeDistribution()
    end
end