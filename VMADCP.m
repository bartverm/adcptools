classdef VMADCP < ADCP &... 
    matlab.mixin.Heterogeneous
% Vessel mounted ADCP wrapper class for adcp structures read with readADCP
%
%   obj=VMADCP() Constructs default object
%
%   obj=VMADCP(...) You can pass different objects to the class on
%   initialization. Depending on the class of the object it will be
%   assigned to a property. On top of those recognized by the ADCP
%   superclass, VMADCP also add support for:
%   - ShipVelocityProvider to the shipvel_provider property
%   
%
%   VMADCP properties (see also inherited properties of ADCP):
%   shipvel_provider - configure how to obtain ship velocity
%
%   VMADCP read-only properties:
%
%   *Data positioning*
%   bt_vertical_range - vertical range to observed bed 
%   slant_range_to_bed - slant range to observed bed
%
%   VMADCP methods
%   bed_offset - xyz offset to the bed
%   bed_position - xyz positions of the bed
%   water_velocity - get velocity profile data corrected for ship motions
%   btvel - ship velocity detected with  bottom tracking
%   plot_bed_position - plot position of detected bed
%
% see also: ADCP

    properties
    % VMADCP/shipvel_provider 
    %
    %   Allows to set how ship velocity is estimated. This is a vector with
    %   ship velocity providers. If ship velocity cannot be obtained from
    %   the first provider, is obtained from the second, and so forth. The
    %   order of the objects in the vector sets the order of preference.
    %   Objects are of class ShipVelocityProvider.
    %   Default is [ShipVelocityFromBT; ShipVelocityFromGPS]
    %
    %   see also: VMADCP, LatLonProvider
        shipvel_provider (:,1) ShipVelocityProvider

    % VMADCP/depth_transducer
    %
    %   Sets the depth of the transducer below the water surface. The value
    %   is always a scalar, non-negative value. The default is 0.
    %
    %   see also: VMADCP
        depth_transducer (1,1) double {mustBeFinite, mustBeReal,...
            mustBeNonnegative} = 0
    end
    properties(Dependent, SetAccess=private)
        % VMADCP/bt_vertical_range read only property.
        %
        % vertical range to the bed detected with bottom tracking.
        %
        % see also: VMADCP
        bt_vertical_range
        
        % VMADCP/slant_range_to_bed
        % 
        % slant range (i.e. along the acoustic beam) to the bed
        %
        % see also: VMADCP
        slant_range_to_bed

    end
    
    methods
        function obj = VMADCP(varargin)
            obj = obj@ADCP(varargin{:})
            obj.water_level_object=ConstantWaterLevel(0);
            obj.vertical_position_provider=ADCPVerticalPositionFromWaterLevel();
            obj.shipvel_provider=[ShipVelocityFromBT; ShipVelocityFromGPS];
            if numel(obj.filters)==1 && isa(obj.filters,'Filter')
                obj.filters=SideLobeFilter;
            end
            for ca=1:nargin
                if isa(varargin{ca},'ShipVelocityProvider')
                    obj.shipvel_provider=varargin{ca};
                end
            end
        end

        function val = get.slant_range_to_bed(obj)
            tm = obj.xform(CoordinateSystem.Instrument, CoordinateSystem.Earth);
            tm(:,:,:,4) = [];
            tm(:,:,4,:) = [];
            beam_m = obj.instrument_matrix_provider.beam_orientation_matrix(obj);
            tm = helpers.matmult(beam_m,tm,3,4);
            tm(:,:,:,1:2)=[];
            val = -obj.bt_vertical_range ./ tm;
        end

        function val = get.bt_vertical_range(obj)
            val = obj.get_bt_vertical_range;
        end
    end

    methods(Sealed)
        function plot_orientations(obj)
            obj.plot_orientations@ADCP;
        end
        function varargout = velocity(obj, varargin)
            varargout = cell(1,min(nargout,numel(obj)));
            [varargout{:}] = obj.velocity@ADCP(varargin{:}); 
        end
        function varargout = xform(obj, varargin)
            varargout = cell(1,min(nargout,numel(obj)));
            [varargout{:}] = obj.xform@ADCP(varargin{:}); 
        end
        function varargout = depth_cell_offset(obj, varargin)
            varargout = cell(1,min(nargout,numel(obj)));
            [varargout{:}] = obj.depth_cell_offset@ADCP(varargin{:}); 
        end        
        function varargout=water_velocity(obj,varargin)
        % VMADCP/water_velocity returns corrected water velocity
        %
        %   vel=water_velocity(obj) returns the velocity corrected for ship
        %   motions in coordinate system specified by obj.coordinate_system
        %
        %   vel=water_velocity(obj,dst) allows to specify a destination
        %   coordinate system
        %
        % see also: VMADCP
            if ~isscalar(obj)
                varargout = obj.run_get_method(nargout,...
                    'water_velocity',varargin{:});
                return
            end
            if nargin < 2
                varargin{1} = obj.coordinate_system;
            end
            varargout{1} = obj.velocity(varargin{1}) +...
                obj.ship_velocity(varargin{1});
        end

        function varargout = btvel(obj, varargin)
        % VMADCP/btvel returns the bottom tracking based ship velocity
        %
        %   btvel=btvel(obj) return the ship velocity as determined with 
        %   the bottom tracking in m/s in the coordinate system specified 
        %   in obj.coordinate_system.
        %
        %   btvel=btvel(obj,dst) allows to specify the destination
        %   coordinate system
        %
        % see also: VMADCP
            if ~isscalar(obj)
                varargout = obj.run_get_method(nargout,...
                    'btvel',varargin{:});
                return
            end
            if nargin < 2
                varargin{1} = obj.coordinate_system;
            end
            varargout{1} = obj.get_btvel(varargin{:});
        end

        function varargout=bed_offset(obj,varargin)
            % Get the offset from the ADCP to the observed bed in m
            %
            %   pos = bed_offset(obj) returns the x,y and z components of the vector
            %   pointing from the ADCP to the detected bed. 
            %
            %   pos = bed_offset(obj, dst) allows to additionaly specify the coordinate
            %   system to get the offsets in as a CoordinateSystem object.
            %
            %   see also: VMADCP, CoordinateSystem, slant_range_to_bed, bed_position
            if ~isscalar(obj)
                varargout = obj.run_get_method(nargout,...
                    'bed_offset', varargin{:});
                return
            end
            if nargin < 2 
                varargin{1}=CoordinateSystem.Earth;
            end
            tm = obj.xform(CoordinateSystem.Instrument,...
                varargin{1}, 'Geometry', true);
            imp = obj.instrument_matrix_provider;
            beam_m = imp.beam_orientation_matrix(obj);
            tm(:,:,:,4)=[];
            tm(:,:,4,:)=[];
            tm = helpers.matmult(beam_m,tm);
            pos=tm.*obj.slant_range_to_bed;
            varargout{1} = pos;
        end

        function varargout=bed_position(obj)
        % Position of the bed in geographic coordinates
        %
        %   pos = bed_position(obj) returns the x,y and z coordinates of the bed as
        %   detected by the bottom tracking of the ADCP.
        %
        % see also: VMADCP, bed_offset
            if ~isscalar(obj)
                varargout = obj.run_get_method(nargout,'bed_position');
                return
            end
            pos = obj.bed_offset(CoordinateSystem.Earth);
            pos(:,:,:,1) = pos(:,:,:,1) + ...
                repmat(obj.horizontal_position(1, :), [1, 1, 4]);
            pos(:,:,:,2) = pos(:,:,:,2) + ...
                repmat(obj.horizontal_position(2, :), [1, 1, 4]);
            pos(:,:,:,3) = pos(:,:,:,3) + ...
                permute(obj.vertical_position, [3, 2, 4, 1]);
            varargout{1} = pos;
        end

        function varargout=ship_velocity(obj,varargin)
            if ~isscalar(obj)
                varargout = obj.run_get_method(nargout,'ship_velocity',...
                    varargin{:});
                return
            end
            varargout{1} = ...
                obj.shipvel_provider.ship_velocity(obj,varargin{:});
        end
 
        function handle_track=plot_track(obj,varargin)
            % plot the track sailed by the vessel
            %
            %   plot_track(obj) plots the track sailed by the vessel
            %
            %   plot_track(obj, ...) pass additional parameters to the plot function 
            %
            %   handle_track = plot_track(...) returns the handle to plot object
            %
            %   see also: VMADCP, plot, plot_all
            for cin=1:numel(varargin)
                if isa(varargin{cin},'matlab.graphics.axis.Axes')
                    ca=varargin{cin};
                    varargin(cin)=[];
                    break
                end
            end
            if ~exist('ca','var')
                ca=gca;
            end
            if ~isscalar(obj)
                hold_stat = get(gca,'NextPlot');
                isef = cellfun(@(x) isa(x,'EnsebleFilter'),varargin);
                if any(isef)
                    warning('Ignoring EnsembleFilter input when plotting an array')
                    varargin(isef)=[];
                end
                for co = 1:numel(obj)
                    obj(co).plot_track(ca,varargin{:})
                end
                set(gca,'NextPlot', hold_stat)
                legend
                return
            end

            arg_used=false(size(varargin));
            ensfilt=EnsembleFilter(obj);
            for cin=1:numel(varargin)
                if isa(varargin{cin},'EnsembleFilter')
                    ensfilt=varargin{cin};
                    arg_used(cin)=true;
                end
            end
            varargin(arg_used) = [];
            ht = nan(numel(ensfilt), 1);
            hold_stat = get(ca, 'NextPlot');
            hold on
            leg = obj.description;
            trk = '';
            for ce = 1 : numel(ensfilt)
                filt =~ ensfilt(ce).all_cells_bad(obj);
                if numel(ensfilt) > 1
                    trk = [' track ', num2str(ce)];
                end
                ht(ce) = plot(ca,...
                    obj.horizontal_position(1,filt),...
                    obj.horizontal_position(2,filt),...
                    varargin{:},'DisplayName',leg + trk);
            end
            if numel(ensfilt) > 1
                legend
            end
            set(gca,'NextPlot',hold_stat)
            axis equal
            xlabel(ca,...
                [obj.horizontal_position_provider.description, ' x (m)']);
            ylabel(ca,...
                [obj.horizontal_position_provider.description, ' y (m)']);
            if nargout > 0
                handle_track = ht;
            end
        end

        function [handle_scat, handle_cbar]=plot_bed_position(obj,varargin)
            pos = cell(size(obj));
            [pos{:}]=obj.bed_position();
            pos = reshape([pos{:}],[], 3);
            hs=scatter3(pos(:,1),pos(:,2),pos(:,3),10,pos(:,3),'filled',varargin{:});
            set(gca,'DataAspectRatio',[1 1 .2])
            view(0,90)
            xlabel([obj(1).horizontal_position_provider.description,' x (m)']);
            ylabel([obj(1).horizontal_position_provider.description,' y (m)']);
            zlabel('z (m)');
            hc=colorbar;
            ylabel(hc,'Bed elevation (m)')
            if nargout > 0
                handle_scat=hs;
            end
            if nargout > 1
                handle_cbar=hc;
            end
        end

        function plot_backscatter(obj)
            hf=plot_backscatter@ADCP(obj);
            add_bed_and_surface(obj,hf,false)
        end

        function plot_track_velocity(obj)
            vel = obj.water_velocity(CoordinateSystem.Earth);
            da_vel = mean(vel,1,"omitnan");
            ph = obj.horizontal_position;
            hold_stat = get(gcf,'NextPlot');
            plot(ph(1,:), ph(2,:),'Color',[.2 .2 .2])
            axis equal
            hold on
            xlabel([obj.horizontal_position_provider.description,' x (m)']);
            ylabel([obj.horizontal_position_provider.description,' y (m)']);
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            xpos = xl(1)+diff(xl)/10;
            ypos = yl(1)+diff(yl)/10;
            text(xpos,ypos,'1  m/s','HorizontalAlignment','left',...
                VerticalAlignment='top',BackgroundColor='w')
            hq=quiver([ph(1,:) xpos],...
                [ph(2,:) ypos], [da_vel(1,:,1) 1], [da_vel(1,:,2) 0]);
            set(gcf,'NextPlot',hold_stat)
        end

        function plot_all(obj)
            plot_all@ADCP(obj)
            figure
            obj.plot_track
            figure
            obj.plot_bed_position
        end

        function add_bed_and_surface(obj,hf,av_beams)
            no = numel(obj);
            if nargin<2
                hf=gcf;
            end
            if nargin < 3
                av_beams=true;
            end
            bed_pos = cell(size(obj));
            [bed_pos{:}] = obj.bed_position;
            bed_pos=cellfun(@(x) squeeze(x(:,:,:,3)), bed_pos, ...
                'UniformOutput', false);
            if av_beams
                bed_pos = cellfun(@(x) repmat(mean(x,2,'omitnan'),...
                    1,size(x,2)), bed_pos, 'UniformOutput',false);
            end
            c=get(hf,'Children');
            if no>1
                assert(isa(c, 'matlab.graphics.layout.TiledChartLayout'),...
                    ['Can only add bed and surface for object arrays',...
                    ' when figure containes tiled layour.']);
                assert(c.GridSize(2) == no, ['For object arrays, ',...
                    'grid size of tiled layout should have as many ',...
                    'columns as number of objects.'])
            end
            if isa(c, 'matlab.graphics.layout.TiledChartLayout')
                c = c.Children;
            end
            cb = size(bed_pos{1},2);
            co = no + 1;
            for cg=c'
                if isa(cg,'matlab.graphics.axis.Axes')
                    co = co - 1;
                    if co < 1
                        co = no;
                        cb = cb - 1;
                        if cb < 1
                            break
                        end
                    end
                    set(cg,'NextPlot','add')
                    plot(cg, obj(co).time, bed_pos{co}(:,cb),...
                        'k','Linewidth',2)
                    plot(cg, obj(co).time, obj(co).water_level,...
                        'b','Linewidth',2)
                    set(cg, 'ylim',...
                        [min(bed_pos{co}(:)) max(obj(co).water_level)])
                end
            end
        end

        function plot_velocity(obj,vel)
            if nargin < 2
                vel = cell(1,numel(obj));
                [vel{:}]=obj.water_velocity(CoordinateSystem.Earth);
            end
            hf=plot_velocity@ADCP(obj,vel);
            add_bed_and_surface(obj,hf)
        end
    end
    methods(Access = protected, Abstract)
        val = get_bt_vertical_range(obj)
        val = get_btvel(obj, dst)
        val = get_velocity(obj,dst)
    end
end