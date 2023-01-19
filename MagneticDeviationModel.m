classdef MagneticDeviationModel < handle
% Provides local magnetic deviation correction
%
%   Abstract class that can be subclassed to define a magnetic deviation
%   model. Subclasses must implement the magnetic_deviation method
%
%   see also: MagneticDeviationConstant, MagneticDeviationTwoCycle
    methods (Abstract)
        % magnetic_deviation returns the magnetic deviation for adcp object
        %
        %   val = magnetic_deviation(adcp) returns the magnetic deviation
        %   in degrees that is added to the heading to obtain the correct
        %   heading
        %
        %   see also: MagneticDeviationModel
        val = magnetic_deviation(adcp)
    end
end

