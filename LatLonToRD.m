classdef LatLonToRD < LatLonToProjection
% Rijksdriehoeks (RD) dutch coordinate system
%
%   LatLonToRD properties:
%   description - 'RD' (Constant property)
%
%   LatLonToRD methods:
%   xy - obtain x,y RD coordinates from lat, lon geographic coordinates
%   ll - obtain lat, lon geographic coordinates from x,y RD coordinates
%
%   see also: VMADCP, LatLonToProjection
    properties
        description = 'RD'
    end
    properties(Access = private, Constant)
    X0      = 155000
    Y0      = 463000
    phi0    = 52.15517440
    lam0    = 5.38720621
    end
    methods
        function obj=LatLonToRD(varargin)
            obj=obj@LatLonToProjection(varargin{:});
        end
        function [lat, lon]=ll(obj, x, y)
            Kp = [0,2,0,2,0,2,1,4,2,4,1];
            Kq = [1,0,2,1,3,2,0,0,3,1,1];
            Kpq = [3235.65389,-32.58297,-0.24750,-0.84978,-0.06550,-0.01709,-0.00738,0.00530,-0.00039,0.00033,-0.00012];

            Lp = [1,1,1,3,1,3,0,3,1,0,2,5];
            Lq = [0,1,2,0,3,1,1,2,4,2,0,0];
            Lpq = [5260.52916,105.94684,2.45656,-0.81885,0.05594,-0.05607,0.01199,-0.00256,0.00128,0.00022,-0.00022,0.00026];

            dX = 1E-5 * ( x - obj.X0 );
            dY = 1E-5 * ( y - obj.Y0 );
        
        lat = 0;
        lon = 0;

        for k = 1:numel(Kpq)
            lat = lat + ( Kpq(k) .* dX.^Kp(k) .* dY.^Kq(k) );
        end
        lat = obj.phi0 + lat / 3600;

        for l = 1:numel(Lpq)
            lon = lon + ( Lpq(l) .* dX.^Lp(l) .* dY.^Lq(l) );
        end
        lon = obj.lam0 + lon / 3600;

        end

        function [x, y]=xy( obj, lat, lon )
        
        Rp = [0,1,2,0,1,3,1,0,2];
        Rq = [1,1,1,3,0,1,3,2,3];
        Rpq = [190094.945,-11832.228,-114.221,-32.391,-0.705,-2.340,-0.608,-0.008,0.148];

        Sp = [1,0,2,1,3,0,2,1,0,1];
        Sq = [0,2,0,2,0,1,2,1,4,4];
        Spq = [309056.544,3638.893,73.077,-157.984,59.788,0.433,-6.439,-0.032,0.092,-0.054];

        dPhi = 0.36 * ( lat - obj.phi0 );
        dLam = 0.36 * ( lon - obj.lam0 );

        x = 0;
        y = 0;

        for r = 1:numel(Rpq)
            x = x + ( Rpq(r) .* dPhi.^Rp(r) .* dLam.^Rq(r) ); 
        end
        x = obj.X0 + x;

        for s = 1:numel(Spq)
            y = y + ( Spq(s) .* dPhi.^Sp(s) .* dLam.^Sq(s) );
        end
        y = obj.Y0 + y;
        end
    end
end