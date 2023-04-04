classdef MagneticDeviationTwoCycle < MagneticDeviationModel
    properties (SetAccess = protected, GetAccess = public)
        % constant offset of compass in degrees
        %
        % see also: MagneticDeviationTwoCycle
        offset = 0


        % sine amplitude of one cycle deviation in degrees
        %
        % see also: MagneticDeviationTwoCycle
        onecycle_sin = 0

        % cosine amplitude of one cycle deviation in degrees
        %
        % see also: MagneticDeviationTwoCycle
        onecycle_cos = 0

        % sine amplitude of two cycle deviation in degrees
        %
        % see also: MagneticDeviationTwoCycle
        twocycle_sin = 0

        % cosine amplitude of two cycle deviation in degrees
        %
        % see also: MagneticDeviationTwoCycle
        twocycle_cos = 0
    end
    properties(Access=protected)
        adcp_head_provider % stores heading provider property of adcp
        % object
    end
    methods
        function val = magnetic_deviation(obj, adcp)
            obj.unset_deviation_correction(adcp)
            head=adcp.heading;
            obj.set_deviation_correction(adcp)

            % compute correction
            val = obj.offset +...
                obj.onecycle_sin * sind(head) + ...
                obj.onecycle_cos * cosd(head) + ...
                obj.twocycle_sin * sind(2 * head) + ...
                obj.twocycle_cos * cosd(2 * head);
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
            d_ang=bt_ang-gps_ang;
            d_ang = atan2d(sind(d_ang), cosd(d_ang));
            head = vmadcp.heading;
            head2 = atan2d(...
                .5*(sind(head(1:end-1))+sind(head(2:end))),...
                .5*(cosd(head(1:end-1))+cosd(head(2:end))));

            obj.set_deviation_correction(vmadcp)

            % fit model
            M = [sind(head2); cosd(head2); sind(2*head2); cosd(2*head2)]';

            pars = robustfit(M,d_ang,'welsch');

            % store model parameters
            obj.offset = pars(1);
            obj.onecycle_sin = pars(2);
            obj.onecycle_cos = pars(3);
            obj.twocycle_sin = pars(4);
            obj.twocycle_cos = pars(5);

            % plot if requested
            if make_plot
                figure
                mod_bias = pars(1) +...
                    pars(2).*M(:,1)'+pars(3).*M(:,2)'+...
                    pars(4).*M(:,3)'+pars(5).*M(:,4)';

                scatter(head2,d_ang,'.')
                hold on
                plot(head2, mod_bias,'r.')
                xlabel('Heading (degrees)')
                ylabel('Heading correction (degrees)')
                legend('Raw differences','Two cycle correction')
            end
        end
    end
    methods(Static)
        function plot_tracks(vmadcp)
            figure
            shipvel = vmadcp.ship_velocity(CoordinateSystem.Earth);
            bt_xvel = shipvel(1,:,1);
            bt_yvel = shipvel(1,:,2);
            bt_xvel2 = .5*(bt_xvel(1:end-1)+bt_xvel(2:end));
            bt_yvel2 = .5*(bt_yvel(1:end-1)+bt_yvel(2:end));
            dt = seconds(diff(vmadcp.time));
            bt_dx = bt_xvel2.*dt;
            bt_dy = bt_yvel2.*dt;
            gps_pos=vmadcp.horizontal_position;
            gps_x = gps_pos(1,:);
            gps_y = gps_pos(2,:);
            bt_x = cumsum([gps_x(1), bt_dx]);
            bt_y = cumsum([gps_y(1), bt_dy]);
            plot(bt_x,bt_y,gps_x,gps_y)
            axis equal
            legend('bt','gps')
            xlabel('x (m)')
            ylabel('y (m)')
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
    end
end