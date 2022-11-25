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
        function vel=water_velocity(obj,varargin)
            vel = water_velocity@VMADCP(obj,varargin{:});

            if nargin > 1
                dst = varargin{1};
            else
                dst = obj.coordinate_system;
            end
            % here we have to deal with bottom track and water track not
            % always using the same set of beams, thus correction in beam
            % coordiantes for ship velocity goes wrong
            src = obj.coordinate_system;
            B = CoordinateSystem.Beam;
            f_corr = src == B & dst == B &...
                obj.raw.WaterTrack.WT_Frequency' ~=...
                obj.raw.BottomTrack.BT_Frequency';

            % we transform from bottom_track beam coordinates, to
            % instrument coordinates and back to water_track coordinates
            bt_b2i = obj.instrument_matrix_provider.b2i_matrix(obj,...
                'BottomTracking',true);
            wt_i2b = obj.instrument_matrix_provider.i2b_matrix(obj);
            tm_bt2wt = helpers.matmult(wt_i2b(:,f_corr,:,:),...
                bt_b2i(:,f_corr,:,:));
            ship_vel = obj.ship_velocity(dst);
            ship_vel(:,f_corr,:,:) = helpers.matmult(tm_bt2wt,...
                ship_vel(:,f_corr,:,:));
            app_vel = obj.velocity(dst);
            vel(:,f_corr,:) = app_vel(:,f_corr,:,:) +...
                ship_vel(:,f_corr,:,:);

        end
    end
    methods(Access = protected)
        function r = get_bt_vertical_range(obj)
            r = permute(obj.raw.BottomTrack.BT_Beam_Depth, [3, 1, 2]);
            r(r==0)=nan;
            r = -r;
            % Untilted sontek ADCP is downlooking, so these should be 
            % negative
        end
        function btvel=get_btvel(obj,dst)
            if nargin < 2
                dst=obj.coordinate_system;
            end
            btvel=permute(obj.raw.BottomTrack.BT_Vel, [3,1,2]);
             if nargin > 1 && ~all(dst == obj.coordinate_system)
                tm=obj.xform(dst,[],'BottomTracking',true);
                btvel=helpers.matmult(tm, btvel);
            end
        end
    end

end