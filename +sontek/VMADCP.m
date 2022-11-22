classdef VMADCP < sontek.ADCP & VMADCP
    methods
        %%% Constructor method %%%
        function obj=VMADCP(varargin)
            obj=obj@sontek.ADCP(varargin{:});
            obj=obj@VMADCP(varargin{:});
            obj.horizontal_position_provider = [
                LatLonToUTM(sontek.LatLonFromGPS)
                ];
        end
    end
    methods(Access = protected)
        function r = get_bt_vertical_range(obj)
            r = permute(obj.raw.BottomTrack.BT_Beam_Depth, [3, 1, 2]);
            r(r==0)=nan;
        end
        function btvel=get_btvel(obj,dst)
            if nargin < 2
                dst=obj.coordinate_system;
            end
            btvel=permute(obj.raw.BottomTrack.BT_Vel, [3,1,2]);
             if nargin > 1 && ~all(dst == obj.coordinate_system)
                tm=obj.xform(dst);
                btvel=helpers.matmult(tm, btvel);
            end
        end
    end

end