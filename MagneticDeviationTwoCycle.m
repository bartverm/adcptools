classdef MagneticDeviationTwoCycle < MagneticDeviationModel
    properties (SetAccess = protected, GetAccess = public)
        pars_c
        pars_s
    end
    properties(Access=protected)
        adcp_head_provider % stores heading provider property of adcp 
        % object
    end
    methods
        function val = magnetic_deviation(obj, adcp)
            obj.unset_deviation_correction(adcp)
            head=adcp.heading;
            obj.set_deviation_correction(adcp);
            val = obj.head_cor(head);
        end

        function estimate_deviation(obj,vmadcp,make_plot)
            if nargin < 3
                make_plot = false;
            end

            obj.unset_deviation_correction(vmadcp)
            
            % compute angles between bottom track and gps tracks
            shipvel = vmadcp.ship_velocity(CoordinateSystem.Earth);
            bt_xvel = shipvel(1,:,1);
            bt_yvel = shipvel(1,:,2);
            bt_xvel2 = .5*(bt_xvel(1:end-1)+bt_xvel(2:end));
            bt_yvel2 = .5*(bt_yvel(1:end-1)+bt_yvel(2:end));
            dt = seconds(diff(vmadcp.time));
            bt_dx = bt_xvel2.*dt;
            bt_dy = bt_yvel2.*dt;
            bt_ang=atan2d(bt_dy,bt_dx);
            gps_pos=vmadcp.horizontal_position;
            gps_x = gps_pos(1,:);
            gps_y = gps_pos(2,:);
            gps_dx = diff(gps_x);
            gps_dy = diff(gps_y);
            gps_ang = atan2d(gps_dy,gps_dx);
            d_ang=-gps_ang+bt_ang;
            d_ang = atan2d(sind(d_ang), cosd(d_ang));
            head = vmadcp.heading;
            head2 = atan2d(...
                .5*(sind(head(1:end-1))+sind(head(2:end))),...
                .5*(cosd(head(1:end-1))+cosd(head(2:end))));
            
            obj.set_deviation_correction(vmadcp)

            % fit model
            M = obj.model_mat(head2);
            obj.pars_c = robustfit(M,cosd(d_ang),'welsch',[],'off');
            obj.pars_s = robustfit(M,sind(d_ang),'welsch',[],'off');

            % plot if requested
            if make_plot
                mod_bias = obj.head_cor(head2);
                figure
                scatter(head2,d_ang,'.')
                hold on
                plot(head2, mod_bias,'r.')
            end
        end
    end
    methods(Access = protected)
        function unset_deviation_correction(obj, adcp)
        % disables previously set MagneticDeviationModels in adcp object

            head_prov=adcp.heading_provider;

             % store heading provider for use in set_deviation_correction
            obj.adcp_head_provider = copy(head_prov);

            % find heading provider in use
            fgood = find(head_prov.has_data(adcp),1,'first');

            % disable the magnetic deviation model
            head_prov(fgood).magnetic_deviation_model =...
                MagneticDeviationConstant(0);
        end
        function set_deviation_correction(obj,adcp)
        % restore orginal heading provider in adcp object
            adcp.heading_provider = obj.adcp_head_provider;
            obj.adcp_head_provider = [];
        end
        function val = head_cor(obj,head)
        % compute correction given heading
            M = obj.model_mat(head);
            mod_bias_c = M * obj.pars_c;
            mod_bias_s = M * obj.pars_s;
            val = reshape(atan2d(mod_bias_s, mod_bias_c),1,[]);
        end
    end
    methods(Access = protected, Static)
        function val = model_mat(head)
            val = [head*0+1;...
                sind(head);...
                cosd(head);...
                sind(2 * head);...
                cosd(2 * head)]';
        end
    end
end