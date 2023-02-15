classdef MagneticDeviationConstant < MagneticDeviationModel
% magnetic deviation is given by a constant value
%
%   obj = MagneticDeviationConstant construct default object
%   obj = MagneticDeviationConstant(val) set the magnetic deviation to val
%       upon construction
%
%   MagenticDeviationConstant properties:
%   mag_deviation - consant magnetic deviation value in degrees
%
% see also: MagneticDeviationModel, MagneticDeviationTwoCycle
    properties
        mag_deviation (1,1) double {mustBeReal, mustBeFinite} = 0;
    end
    methods
        function obj = MagneticDeviationConstant(val)
            obj = obj@MagneticDeviationModel;
            if nargin > 0
                obj.mag_deviation = val;
            end
        end
        function val = magnetic_deviation(obj, ~)
            val = obj.mag_deviation;
        end
    end
end