classdef VMADCP < ADCP
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

    end
    properties (SetObservable)
        water_level_object (1,1) WaterLevel= ConstantWaterLevel(0);
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

        water_level
    end
    methods
        %%% Constructor method %%%
        function obj=VMADCP(varargin)
            obj=obj@ADCP(varargin{:});
            obj.water_level_object=ConstantWaterLevel(0);
            obj.vertical_position_provider=ADCPVerticalPositionFromWaterLevel(obj.water_level_object, obj);
            obj.horizontal_position_provider=[ProjectedCoordinatesFromViseaExtern; LatLonToUTM];
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
        
        %%% Set and Get methods %%%
        function r=get.bt_vertical_range(obj)
            r=shiftdim(double(obj.raw.btrange)/100,-1);
            r(r==0)=nan;
        end
        function r=get.slant_range_to_bed(obj)
            bangle=obj.beam_angle;
            r=obj.bt_vertical_range./cosd(bangle);
        end
        
        %%% Ordinary methods %%%
        function vel=ship_velocity(obj,dst)
            vel=obj.shipvel_provider.ship_velocity(obj,dst);
        end
        
        function vel=water_velocity(obj,dst)
        % VMADCP/water_velocity returns corrected water velocity
        %
        %   vel=water_velocity(obj) returns the velocity corrected for ship
        %   motions in coordinate system specified by obj.coordinate_system
        %
        %   vel=water_velocity(obj,dst) allows to specify a destination
        %   coordinate system
        %
        % see also: VMADCP
            if nargin < 2
                dst=obj.coordinate_system;
            end
            vel=obj.velocity(dst)+obj.ship_velocity(dst);
        end
        
        function btvel=btvel(obj,dst)
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
            if nargin < 2
                dst=obj.coordinate_system;
            end
            btvel=shiftdim(VMADCP.int16_to_double(obj.raw.btvel)/1000,-1);
             if nargin > 1 && ~all(dst == obj.coordinate_system)
                tm=obj.xform(dst);
                btvel=helpers.matmult(tm, btvel);
            end
        end
        
        function pos=bed_offset(obj,dst)
            % Get the offset from the ADCP to the observed bed in m
            %
            %   pos = bed_offset(obj) returns the x,y and z components of the vector
            %   pointing from the ADCP to the detected bed. 
            %
            %   pos = bed_offset(obj, dst) allows to additionaly specify the coordinate
            %   system to get the offsets in as a CoordinateSystem object.
            %
            %   see also: VMADCP, CoordinateSystem, slant_range_to_bed, bed_position
            if nargin < 2 
                dst=CoordinateSystem.Earth;
            end
            tm=-obj.xform(CoordinateSystem.Beam, dst,'UseTilts',true);
            tm(:,:,:,4)=[];
            pos=tm.*obj.slant_range_to_bed;
        end
        function pos=bed_position(obj)
        % Position of the bed in geographic coordinates
        %
        %   pos = bed_position(obj) returns the x,y and z coordinates of the bed as
        %   detected by the bottom tracking of the ADCP.
        %
        % see also: VMADCP, bed_offset
            pos=obj.bed_offset(CoordinateSystem.Earth);
            pos(:,:,:,1)=pos(:,:,:,1)+repmat(obj.horizontal_position(1,:),[1,1,4]);
            pos(:,:,:,2)=pos(:,:,:,2)+repmat(obj.horizontal_position(2,:),[1,1,4]);
            pos(:,:,:,3)=pos(:,:,:,3)+permute(obj.vertical_position,[3,2,4,1]);
        end
        function wl=get.water_level(obj)
            wl=obj.water_level_object.get_water_level(obj.time);
        end
        function pos=depth_cell_position(obj)
            % Get the geographic position of the depth cells
            %
            %   pos = depth_cell_position(obj) returns the x, y and z coordinates of
            %   the depth cells of the ADCP.
            %
            %   see also: VMADCP, depth_cell_offset
            pos=obj.depth_cell_offset(CoordinateSystem.Earth);
            pos(:,:,:,1)=pos(:,:,:,1)+repmat(obj.horizontal_position(1,:),[1,1,4]);
            pos(:,:,:,2)=pos(:,:,:,2)+repmat(obj.horizontal_position(2,:),[1,1,4]);
            pos(:,:,:,3)=pos(:,:,:,3)+permute(obj.vertical_position,[3,2,4,1]);
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
                if isa(varargin{cin},'matlab.graphics.axis.Axes')% || isa(obj.axes,'matlab.ui.control.UIAxes')
                    ca=varargin{cin};
                    varargin(cin)=[];
                    break
                end
            end
            if ~exist('ca','var')
                ca=gca;
            end
            arg_used=false(size(varargin));
            ensfilt=EnsembleFilter(obj);
            for cin=1:numel(varargin)
                if isa(varargin{cin},'EnsembleFilter')
                    ensfilt=varargin{cin};
                    arg_used(cin)=true;
                end
            end
            varargin(arg_used)=[];
            ht=nan(numel(ensfilt),1);
            hold_stat=get(ca,'NextPlot');
            hold on
            for ce=1:numel(ensfilt)
                filt=~ensfilt(ce).all_cells_bad(obj);
                ht(ce)=plot(ca,obj.horizontal_position(1,filt),obj.horizontal_position(2,filt),varargin{:});
            end
            set(gca,'NextPlot',hold_stat)
            axis equal
            xlabel(ca,[obj.horizontal_position_provider.description,' x (m)']);
            ylabel(ca,[obj.horizontal_position_provider.description,' y (m)']);
            if nargout > 0
                handle_track=ht;
            end
        end
        function [handle_scat, handle_cbar]=plot_bed_position(obj,varargin)
            pos=obj.bed_position();
            pos=reshape(pos,[],3);
            hs=scatter3(pos(:,1),pos(:,2),pos(:,3),10,pos(:,3),'filled',varargin{:});
            set(gca,'DataAspectRatio',[1 1 .2])
            view(0,90)
            xlabel([obj.horizontal_position_provider.description,' x (m)']);
            ylabel([obj.horizontal_position_provider.description,' y (m)']);
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
        function plot_velocity(obj,vel)
            if nargin < 2
                vel=obj.water_velocity(CoordinateSystem.Earth);
            end
            hf=plot_velocity@ADCP(obj,vel);
            add_bed_and_surface(obj,hf)
        end
        function plot_backscatter(obj)
            hf=plot_backscatter@ADCP(obj);
            add_bed_and_surface(obj,hf,false)
        end
        function plot_all(obj)
            plot_all@ADCP(obj)
            figure
            obj.plot_track
            figure
            obj.plot_bed_position
        end
        function add_bed_and_surface(obj,hf,av_beams)
            if nargin<2
                hf=gcf;
            end
            if nargin < 3
                av_beams=true;
            end
            bed_pos=obj.bed_position;
            bed_pos=squeeze(bed_pos(:,:,:,3));
            if av_beams
                bed_pos=repmat(mean(bed_pos,2,'omitnan'), 1,size(bed_pos,2));
            end
            c=get(hf,'Children');
            cb=size(bed_pos,2)+1;
            for cg=c'
                if isa(cg,'matlab.graphics.axis.Axes')
                    cb=cb-1;
                    set(cg,'NextPlot','add')
                    t=obj.time;
                    t=seconds(t-t(1));
                    plot(cg,t,bed_pos(:,cb),'k','Linewidth',2)
                    set(cg,'ylim',[min(bed_pos(:)) 0])
                end
            end
        end
    end
end