classdef PhasedArrayTransducer < handle
    properties
        frequency(1,1) double {mustBeFinite, mustBePositive} 
        radius(1,1) double {mustBeFinite, mustBePositive}
    end
end