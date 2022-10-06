classdef VaryingWaterLevel < WaterLevel
% Implements a water level varying in time
%
%   This object uses a time series of water level to provide water levels
%   at a given time by linear interpolation. It does not extrapolate.
%
%   obj=VaryingWaterLevel() constructs and empty water level object
%   obj=VaryingWaterLevel(...) pass arguments to constructor that are
%   assigned to properties based on their class as follows:
%   - datetime object are assigned to the time property
%   - double arrays are passed to the water_level property
%
%   VaryingWaterLevel properties:
%   time - time of the water level measurements as datetime array
%   water_level - water level as double array
%
%   VaryingWaterLevel methods:
%   plot - plots the water level time series
%
%   see also: VMADCP, ADCPVerticalPositionFromWaterLevel, WaterLevel
    properties
% VaryingWaterLevel/time time of water level measurements
%
%   time is a datetime row vector holding the time the water level
%   measurements were made. It should have the same number of elements as
%   the water_level property for the class to function correctly
%
%   see also: VaryingWaterLevel, water_level, datetime
        time (1,:) datetime = datetime()
        constituents(1,:) cell = {};
        names(1,:) cell = {};
        parameters(1,:) double = []
    
% VaryingWaterLeve/water_level water level series
%
%   water_level is a double row vector holding the measured water levels in
%   m. It should have the same number of elements as the time property for
%   the class to function properly.
%
%   see also: VaryingWaterLevel, time
        water_level (1,:) double {mustBeReal, mustBeFinite} = 0
    end
    methods
        function obj=VaryingWaterLevel(varargin)
            for ca=1:nargin
                if isa(varargin{ca},'double')
                    obj.water_level=varargin{ca};
                elseif isa(varargin{ca},'datetime')
                    obj.time = varargin{ca};
                else
                    warning(['Unused input number ',num2str(ca),' of type ',class(varargin{ca})])
                end
            end
        end
        function val=get_water_level(obj,time)
            val=interp1(obj.time,obj.water_level,time,'linear');
        end

        function M = get_model(obj)
            % Fitting tidal constituents to waterlevel
            npars = 2*numel(obj.constituents) + 1;
            d_hours = (datenum(obj.time)-datenum(obj.time(1)))*24;
            T = const_to_periods(obj.constituents);
            M = ones(numel(d_hours), npars);
            for c = 1:numel(obj.constituents)
                M(:,2*c) = cos(2*pi/T(c) * d_hours);
                M(:,2*c + 1) = sin(2*pi/T(c) * d_hours);
            end
        end

        function obj = get_parameters(obj)
            M = obj.get_model();
            obj.parameters = (M'*M)\M'*obj.water_level';
            obj.names = obj.get_names();
        end

       function names = get_names(obj)
            names = cell(1,2*numel(obj.constituents)+1);
            names{1} = 'eta: M0A';
            for i = 1:numel(obj.constituents)
                names{2*i} = sprintf('%s%s%s', 'eta: ', obj.constituents{i},'A');
                names{2*i+1} = sprintf('%s%s%s', 'eta: ', obj.constituents{i},'B');
            end
            obj.names = names;
        end

        function hp=plot(obj,varargin)
% plots the water level time series
%
%   obj.plot() plots the water level time series
%   obj.plot(...) pass additional arguments to the plot function, e.g. 'r.'
%   to plot red dots etc.
%   hp = obj.plot() return the handle to the water level time series
%
%   see also: VaryingWaterLevel
            htmp=plot(obj.time,obj.water_level,varargin{:});
            xlabel('Time')
            ylabel('Water level (m)')
            if nargout>0
                hp=htmp;
            end
        end
    end
end
