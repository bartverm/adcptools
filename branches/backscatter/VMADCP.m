classdef VMADCP < ADCP
% Vessel mounted ADCP wrapper class for adcp structures read with readADCP
%
%   obj=VMADCP() Constructs default object
%
%   obj=VMADCP(raw) Allows to assign the adcp structure to be used to the 
%   raw property
%
%   VMADCP properties:
%   raw - adcp structure read by readADCP.m
%   filters - filters for profiled data. Defaults to SideLobeFilter
%   timezone - timezone of the data
%   type - type of ADCP being used
%   xy_cor_system - geographic coordinate system to be used
%   adcp_elevation - elevation of the ADCP
%
%   VMADCP read-only properties:
%   *Instrument characteristics*
%   transducer - object describing transducer characteristics
%   is_workhorse - return whether ADCP is a Workhorse ADCP
%
%   *Sizing*
%   fileid - ID of file ensemble was read from
%   time - time ensembles were measured
%   nensembles - number of ensembles
%   nbeams - number of beams
%   ncells - number of cells
%
%   *Coordinate systems properties*
%   coordinate_system - coordinate system used
%   beam_angle - angle of acoustic beam with vertical
%   convexity - convexity of the instrument
%   three_beam_solutions_used - whether three beams solutions were used
%   tilts_used_in_transform - whether tilts were used in transformations
%   bin_mapping_used - wheteher bin mapping was used during measurements     
%
%   *Orientation*
%   pitch - pitch rotation angle
%   roll - roll rotation
%   heading - heading rotation
%   headalign - heading alignment 
%   is_upward - whether instrument was pointing upward
%
%   *Data positioning*
%   cellsize - vertical size of depth cell
%   lengthxmitpulse - length of transmitted pulse
%   blanking - blanking distance
%   distmidfirstcell - vertical distance to the first measured depth cell
%   depth_cell_slant_range - slant range to depth cell
%   bt_vertical_range - vertical range to observed bed 
%   slant_range_to_bed - slant range to observed bed
%
%   *Backscatter*
%   current - transmit current of transducer (A)
%   current_factor - factor for current computation from ADC channel
%   voltage - transmit voltage of transducer (V)
%   voltage_factor - factor for voltage computation from ADC channel
%   power - transmit power of transducer (W)
%   attitude_temperature - temperature of transducer (Celsius)
%   intensity_scale - intensity scale factor (dB/m)
%   echo - Raw echo intennstiy (dB)
%   backscatter - Volumne backscatter strength (dB)
%
%   VMADCP methods
%   bad - filter in use                           
%   xform - coordinate transformation matrices
%   depth_cell_offset - xyz offset from transducer to depth cell
%   depth_cell_position - xyz positions of depth cells in projected coords.
%   bed_offset - xyz offset to the bed
%   bed_position - xyz positions of the bed
%   lat_lon - lat lon positions of ADCP
%   projected - xy projected coordinates of the ADCP
%   velocity - get velocity profile data
%   water_velocity - get velocity profile data corrected for ship motions
%   btvel - ship velocity detected with  bottom tracking
%   plot_filters - plot active profile data filters
%   plot_orientations - plot orientations of ADCP
%   plot_velocity - plot velocity profiling data
%   plot_track - plot track of the ADCP
%   plot_bed_position - plot position of detected bed
%
%   VMADCP static methods:
%   beam_2_instrument - beam to instrument transformation matrices
%   head_tilt - tilt matrices 
%   instrument_2_beam - instrument to beam transformation matrices
%   int16_to_double - transform int16 to double values
%   invert_xform - invert transformation matrices
%
% see also: ADCP
    properties
    % VMADCP/xy_cor_system 
    %
    %   Projected geographic coordinate system to be used. Default is the
    %   UTMCoordinateSystem. Value must be an object of a class inherited
    %   from ProjectedCoordinateSystem.
    %
    %   see also: VMADCP
        xy_cor_system (1,1) ProjectedCoordinateSystem = UTMCoordinateSystem;
    
    % VMADCP/adcp_elevation
    %
    %   specifies the elevation of the ADCP during the measurement. This
    %   can be used to correct for water level changes during data
    %   collection
        adcp_elevation (1,:) double {mustBeFinite} = 0;
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
        %%% Constructor method %%%
        function obj=VMADCP(varargin)
            obj=obj@ADCP(varargin{:});
            obj.filters=SideLobeFilter;
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
            vel=obj.velocity(dst)-obj.btvel(dst);
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
        function [lat, lon]=lat_lon(obj,ll_provider)
            if nargin<2
                ll_provider=LatLonVisea;
            end
            [lat,lon]=ll_provider.lat_lon(obj.raw);
        end
        function [x, y]=projected(obj)
            [lat, lon]=obj.lat_lon;
            [x, y]= obj.xy_cor_system.xy(lat,lon);
        end
        function pos=bed_offset(obj,dst)
            if nargin < 2 
                dst=CoordinateSystem.Earth;
            end
            tm=-obj.xform(CoordinateSystem.Beam, dst);
            tm(:,:,:,4)=[];
            pos=tm.*obj.slant_range_to_bed;
        end
        function pos=bed_position(obj)
            pos=obj.bed_offset(CoordinateSystem.Earth);
            [x,y]=obj.projected();
            pos(:,:,:,1)=pos(:,:,:,1)+x;
            pos(:,:,:,2)=pos(:,:,:,2)+y;
            pos(:,:,:,3)=pos(:,:,:,3)+obj.adcp_elevation;
        end
        function pos=depth_cell_position(obj)
            pos=obj.depth_cell_offset(CoordinateSystem.Earth);
            [x,y]=obj.projected();
            pos(:,:,:,1)=pos(:,:,:,1)+x;
            pos(:,:,:,2)=pos(:,:,:,2)+y;
            pos(:,:,:,3)=pos(:,:,:,3)+obj.adcp_elevation;
        end
        function plot_track(obj)
            [x,y]=obj.projected;
            plot(x,y,'.-');
            axis equal
            xlabel([obj.xy_cor_system.description,' x (m)']);
            ylabel([obj.xy_cor_system.description,' y (m)']);
        end
        function plot_bed_position(obj)
            pos=obj.bed_position();
            pos=reshape(pos,[],3);
            scatter3(pos(:,1),pos(:,2),pos(:,3),10,pos(:,3),'filled');
            set(gca,'DataAspectRatio',[1 1 .2])
            view(0,90)
            xlabel([obj.xy_cor_system.description,' x (m)']);
            ylabel([obj.xy_cor_system.description,' y (m)']);
            zlabel('z (m)');
            hc=colorbar;
            ylabel(hc,'Bed elevation (m)')
        end
        function plot_velocity(obj,vel)
            if nargin < 2
                vel=obj.water_velocity(CoordinateSystem.Earth);
            end
            plot_velocity@ADCP(obj,vel)
            add_bed_and_surface(obj,gcf)
        end
        function plot_backscatter(obj)
            plot_backscatter@ADCP(obj)
            add_bed_and_surface(obj,gcf)
        end
        function add_bed_and_surface(obj,hf)
            if nargin<2
                hf=gcf;
            end
            bed_pos=obj.bed_offset;
            bed_pos=squeeze(bed_pos(:,:,:,3));
            bed_pos=nanmean(bed_pos,2);
            c=get(hf,'Children');
            for cg=c'
                if isa(cg,'matlab.graphics.axis.Axes')
                    set(cg,'NextPlot','add')
                    t=obj.time;
                    t=seconds(t-t(1));
                    plot(cg,t,bed_pos,'k','Linewidth',2)
                    set(cg,'ylim',[min(bed_pos) 0])
                end
            end
        end
    end
end