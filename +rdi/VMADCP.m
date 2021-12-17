classdef VMADCP < VMADCP & rdi.ADCP

    methods
        %%% Constructor method %%%
        function obj=VMADCP(varargin)
            obj=obj@rdi.ADCP(varargin{:});
            obj=obj@VMADCP(varargin{:});
            obj.horizontal_position_provider=[ProjectedCoordinatesFromViseaExtern; LatLonToUTM];
        end
    end
    methods(Access = protected)
        function r=get_bt_vertical_range(obj)
            r=shiftdim(double(obj.raw.btrange)/100,-1);
            r(r==0)=nan;
        end
        function btvel=get_btvel(obj,dst)
            if nargin < 2
                dst=obj.coordinate_system;
            end
            btvel=shiftdim(obj.int16_to_double(obj.raw.btvel)/1000,-1);
             if nargin > 1 && ~all(dst == obj.coordinate_system)
                tm=obj.xform(dst);
                btvel=helpers.matmult(tm, btvel);
            end
        end
    end
end