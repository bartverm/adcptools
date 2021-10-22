classdef ConstantWaterLevel < WaterLevel
% Defines a constant water level
%
%   obj = ConstantWaterLevel() constructs a default water level object with
%   the level set to 0
%
%   obj = ConstantWaterLevel(level) constructs the object with the given
%   level
%
%   ConstantWaterLevel properties:
%   level - the water level in m
%
%   ConstantWaterlevel methods:
%   get_water_level - returns the water level
%
%   see also: VMADCP
    properties
% ConstantWaterLevel/level property
%
%   The water level in m. This is a scalar double. Default is 0.
%
%   see also: ConstantWaterLevel
        level(1,1) double {mustBeFinite, mustBeReal} = 0;
    end
    methods
        function obj=ConstantWaterLevel(varargin)
            if nargin > 0
                obj.level=varargin{1};
            end
        end
        function wl=get_water_level(obj,time)
% Return the water level at given time
%
%   wl = get_water_level() returns the water level
%
%   wl = get_water_level(time) returns an array with the same size as time
%   holding the water level
%
%   see also: ConstantWaterLevel
            wl=ones(size(time))*obj.level;
        end
    end
end