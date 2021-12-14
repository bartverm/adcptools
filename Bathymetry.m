classdef Bathymetry < handle
% Abstract class defining how bathymetry is provided to ADCPtools
%
%   Subclasses need to implement the 'get_bed_level' method
%
%   Bathymetry construction (only callable on concrete subclasses):
%   obj=Bathymetry() creates an emtpy bathymetry object
%   obj=Bathymetry(wl) any argument passed of class WaterLevel is assigned
%       to the water_level property. If more than one is passed, the last
%       one is used.
%
%   Bathymetry properties:
%   water_level - defines the water level
%
%   Bathymetry methods:
%   get_depth - compute depth at a given position and time
%   plot - calls the plot function on all concrete array elements 
%
%   See also: BathymetryScatteredPoints
    properties
% Bathymetry/water_level
%
%   scalar object of class WaterLevel, defining the water level. Needed to
%   compute depth given position and time
%
%   see also: Bathymetry, get_depth, WaterLevel
        water_level (1,1) WaterLevel=ConstantWaterLevel;
    end
    methods (Abstract)
% Bathymetry/get_bed_elev abstract method
%
%   z=obj.get_bed_elev(x,y) returns the bed elevation at the locations
%   given in x,y. z has dimensions equal to dimensions of x and y.
%
%   method needs to be reimplemented in concrete subclass
%
%   see also: Bathymetry, get_depth
        z=get_bed_elev(obj,x,y)
    end
    methods
        function obj=Bathymetry(varargin)
            for count_var=1:nargin
                cvar=varargin{count_var};
                if isa(cvar,'WaterLevel')
                    obj.water_level=cvar;
                end
            end
        end


        function d=get_depth(obj,x,y,time)
% Bathymetry/get_depth
%
%   d=obj.get_depth(x,y,time) compute the depth given the position x,y and
%   time. d has the same dimensions as the dimensions of x,y and time.
%
%   see also: Bathymetry, get_bed_elev
            if nargin > 3
                d=obj.water_level.get_depth(obj.get_bed_elev(x,y),time);
            else
                d=obj.obj.water_level.get_depth(obj.get_bed_elev(x,y));
            end
        end

        function plot(obj,varargin)
% Bathymetry/plot
%
%   obj.plot(...) given an array of objects obj, calls the plot function on
%       all concrete array elements in obj making sure all bathymetries are
%       added to same axes. plot method in subclasses should call this plot
%       method if obj is an array, and plot if obj is scalar.
%
%   see also: Bathymetry
            obj(1).plot;
            hold_stat=get(gca,'NextPlot');
            hold on;
            if ~isscalar(obj)
                for ce=2:numel(obj)
                    obj(ce).plot(varargin{:})
                end % for
            end % if
            set(gca,'NextPlot',hold_stat);
        end % function
    end % methods
end % classdef