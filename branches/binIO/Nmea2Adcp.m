classdef Nmea2Adcp < handle
    properties
        pd0
        nmea
    end
    methods
        function obj=Nmea2Adcp(varargin)
            narginchk(0,2);
            if nargin>0
                obj.pd0=varargin{1};
            end
            if nargin>1
                obj.nmea=varargin{2};
            end
        end
        function set.pd0(obj,val)
            validateattributes(val,{'PD0'},{'scalar'},'set.pd0','pd0',2);
            obj.pd0=val;
        end
        function set.nmea(obj,val)
            validateattributes(val,{'NMEA'},{'scalar'},'set.nmea','nmea',2);
            obj.nmea=val;
        end
        
    end
    methods (Abstract)
        data=match(obj)
    end
end