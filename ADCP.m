classdef ADCP < handle
% Abstract base class to represent ADCP datasets
%
%   A = ADCP(...) based on the class of the passed arguments
%   configuration properties are set:
%   Filter - filter
%   HorizontalPositionProvider - horizontal_position_provider
%   VerticalPositionProvider - vertical_position_provider
%
%   ADCP properties (configuration):
%   filters - filters to exclude data from profiles
%   timezone - specifies the timezone used
%   horizontal_position_provider - provides horizontal positioning
%   vertical_position_provider - provides vertical positioning
%   tilts_provider - provides pitch and roll angles
%   heading_provider - provides heading angle
%   instrument_matrix_provider - provides instrument to beam matrix
%   transducer - transducer acoustic properties
%   water - water acoustic properties
%
%   ADCP data properties (read only):
%   nbeams - number of acoustic beams
%   nensembles - number of ensembles
%   ncells - number of depth cells
%   coordinate_system - coordinate system of velocity data
%   beam_angle - angle of acoustic beams with vertical (degrees)
%   pitch - pitch angle (degrees)
%   roll - roll angle (degrees)
%   heading - heading angle (degrees)
%   cellsize - depth cell size (m)
%   time - time and date of ensembles
%   distmidfirstcell - distance to center of first depth cell (m)
%   depth_cell_slant_range - distance along ac. beam to depth cells (m)
%   temperature - instrument temperature (Celsius)
%   salinity - salinity of water (psu)
%   pressure - pressure of water (Pa)
%   echo - received echo intensity (dB)
%   backscatter - volume backscatter strength (dB)
%   horizontal_position - horizontal position of the instrument (m)
%   vertical_position - vertical position of the instrument (m)
%   beam_2_instrument_matrix - beam to instrument transformation matrix
%   instrument_2_beam_matrix - instrument to beam transformation matrix
%
%   ADCP methods:
%   bad - return a mask with ones for data marked as being bad
%   plot_orientations - plots pitch, roll and heading of the ADCP
%   plot_velocity - plot the velocity in Earth coordinates
%   plot_backscatter - plot the volume backscatter strength
%   plot_all - plot all available plots
%   depth_cell_offset - get east,north,up offset of depth cells from ADCP
%   depth_cell_position - get east,north,up position of depth cells
%   
%   ADCP methods (abstract):
%   xform - get transformation matrices
%   velocity - get velocity

    properties
        % ADCP/filters
        %
        % Filters to exclude data from profile data. The filters are given
        % as a vector of object that are subclasses of Filter
        %
        % see also: ADCP, Filter
        filters(:,1) Filter {mustBeNonempty}=Filter;
        
        % ADCP/timezone
        %
        %   Specify the timezone of the data as char. Default is '', i.e.
        %   unzoned. Examples are 'UTC', 'Europe/Amsterdam'
        %
        % see also: ADCP, datetime/TimeZone, timezones
        timezone(1,:) char = ''

        % ADCP/horizontal_position_provider
        %
        %   ADCPHorizontalPosition object specifying the position of the ADCP.
        %
        % see also: ADCP
        horizontal_position_provider(1,:) ADCPHorizontalPosition = ADCPFixedHorizontalPosition;
        
        % ADCP/vertical_position_provider
        %
        %   ADCPVerticalPosition object specifying the position of the ADCP.
        %
        % see also: ADCP
        vertical_position_provider(1,1) ADCPVerticalPosition = ADCPFixedVerticalPosition;

        % ADCP/tilts_provider
        %
        % TiltsProvider object which returns the heading and tilts of the 
        % ADCP.
        %
        % see also: ADCP
        tilts_provider(:,1) TiltsProvider = rdi.TiltsInternal;

        % ADCP/heading_provider
        %
        % HeadingProvider object which returns the heading and tilts of the 
        % ADCP.
        %
        % see also: ADCP
        heading_provider(:,1) HeadingProvider = rdi.HeadingInternal


        % ADCP/instrument_matrix_provider
        %
        %   Specifies the sources for the transformation matrix as a
        %   InstrumentMatrixProvider. This is a vector which allows to
        %   define different sources. The first object capable of providing
        %   the matrix will be used.
        %
        % see also: ADCP, InstrumentMatrixProvider
        instrument_matrix_provider (:,1) InstrumentMatrixProvider =...
            rdi.InstrumentMatrixFromBAngle;
    end
    properties(Access=protected)
        override_transducer
        override_water
    end
    properties(Dependent)
        % ADCP/tranducer
        %
        % returns an acoustics.PistonTransducer object. This object can be
        % modified, and will be used in computations. To reset object to
        % its initial values use the reset_transducer method. Setting the
        % raw property, or the water property will reset the tranducer
        % property
        %
        % see also: ADCP, acoustics.PistonTransducer
        transducer(:,1) acoustics.PistonTransducer 
        
        % ADCP/water
        %
        %   acoustics.Water object specifying the Water characteristics.
        %   use reset_water method to reset the water property to the
        %   values read in the raw adcp data. Changing the raw property
        %   will also reset the water property.
        %
        % see also: ADCP, acoustics.Water
        water(1,1) acoustics.Water;    
    end
    properties(Dependent, SetAccess = private)
        % Number of beams in use by the instrument.
        %
        % see also: ADCP
        nbeams
        
        % Number of ensembles.
        %
        % see also: ADCP
        nensembles
        
        % Number of depth cells.
        %
        % see also: ADCP
        ncells
        
        % Coordinate system used. Returns CoordinateSystems objects.
        %
        % see also: ADCP, CoordinateSystem
        coordinate_system
        
        % Angle of acoustic beam makes with vertical in degrees.
        %
        % see also: ADCP
        beam_angle
        
        % Pitch angle in degrees.
        %
        % see also: ADCP
        pitch
        
        % Roll angle in degrees.
        %
        % see also: ADCP
        roll
        
        % Heading angle in degrees.
        %
        % see also: ADCP
        heading

        % Vertical size of the depth cells (m).
        %
        % see also: ADCP
        cellsize

        % Time the ensembles were measured. Returns datetime objects.
        %
        % see also: ADCP, datetime
        time
        
        % Vertical distance to the middle of the first cell (m).
        %
        % see also: ADCP
        distmidfirstcell
        
        % Slant range, i.e. distance along acoustic beam to depth cells
        %
        % see also: ADCP
        depth_cell_slant_range
        
        % Temperature of the instrument (Celsius)
        %
        % see also: ADCP
        temperature
        
        %  salinity (psu)
        %
        % see also: ADCP
        salinity
        
        % pressure (Pa)
        %
        % see also: ADCP
        pressure

        % Received raw echo intensity (dB).
        %
        % see also: ADCP
        echo
        
        % Received Volume Backscatter strength (dB) computed according to
        % Deines (1999) with the corrections of Gostiaux and van Haren and 
        % the correction in the FSA-031. The backscatter strength is not 
        % corrected for attenuation due to sediment.
        %
        % see also: ADCP, acoustics, Sv2SSC
        backscatter
        
        %   2xN matrix holding the x and y coordinates of the
        %   ADCP in m.
        %
        %   see also: adcp, vertical_position, horizontal_position_provider
        horizontal_position
        
        %   1xN vector holding the z coordinates of the ADCP in m.
        %
        %   see also: adcp, horizontal_position, vertical_position_provider
        vertical_position

        % matrix transforming beam coordinates to instrument coordinates
        %
        %   see also: ADCP, instrument_matrix_provider
        beam_2_instrument_matrix

        % matrix transforming instrument coordinates to beam coordinates
        %
        %   see also: ADCP, instrument_matrix_provider
        instrument_2_beam_matrix

        beam_orientation_matrix
    end
    methods
        %%% Constuctor
        function obj = ADCP(varargin)
            obj.filters=Filter;
            for ca=1:nargin            
                if isa(varargin{ca},'Filter')
                    obj.filters=[obj.filters; varargin{ca}];
                elseif isa(varargin{ca},'ADCPVerticalPosition')
                    obj.vertical_position_provider=varargin{ca};
                elseif isa(varargin{ca},'ADCPHorizontalPosition')
                    obj.horizontal_position_provider=varargin{ca};
                end
            end
        end
        
        %%% Set and get methods
        function val = get.water(obj)
            if isempty(obj.override_water)
                val=acoustics.Water;
                val.temperature=obj.temperature;
                val.salinity=obj.salinity*1000;
            else
                val=obj.override_water;
            end
        end
        function set.water(obj, val)
            if isempty(val)
                obj.water = acoustics.Water.empty;
                return
            end
            assert(isa(val, 'acoustics.Water'), 'Must be an acoustics.Water object')
            assert(isscalar(val), 'Value must be scalar')
            obj.override_water = val;
        end
        function val = get.transducer(obj)
            if isempty(obj.override_transducer)
                val = obj.get_transducer;
            else
                val = obj.override_transducer;
            end
        end
        function set.transducer(obj, val)
            if isempty(val)
                obj.transducer = acoustics.PistonTransducer.empty;
                return
            end
            assert(isa(val, 'acoustics.Transducer'), 'Must be an acoustics.Transducer object')
            assert(isscalar(val), 'Value must be scalar')
            obj.override_transducer = val;
        end
        function val = get.nbeams(obj)
            val = obj.get_nbeams;
        end
        function val = get.nensembles(obj)
            val = obj.get_nensembles;
        end
        function val = get.ncells(obj)
            val = obj.get_ncells;
        end
        function val = get.coordinate_system(obj)
            val = obj.get_coordinate_system;
        end
        function val = get.beam_angle(obj)
            val = obj.get_beam_angle;
        end
        function val = get.beam_orientation_matrix(obj)
            val = obj.instrument_matrix_provider.beam_orientation_matrix(obj);
        end
        function val = get.pitch(obj)
            val = obj.tilts_provider.pitch(obj);
        end
        function val = get.roll(obj)
            val = obj.tilts_provider.roll(obj);
        end
        function val = get.heading(obj)
            val = obj.heading_provider.heading(obj);
        end
        function val = get.cellsize(obj)
            val = obj.get_cellsize;
        end
        function val = get.time(obj)
            val = obj.get_time;
        end
        function val = get.distmidfirstcell(obj)
            val = obj.get_distmidfirstcell;
        end
        function rng=get.depth_cell_slant_range(obj)
            bangle=obj.beam_angle;
            rng=(obj.distmidfirstcell+reshape(0:max(obj.ncells)-1,[],1).*obj.cellsize)./cosd(bangle);
        end
        function val = get.temperature(obj)
            val = obj.get_temperature;
        end
        function val = get.salinity(obj)
            val = obj.get_salinity;
        end
        function val = get.pressure(obj)
            val = obj.get_pressure;
        end
        function val = get.echo(obj)
            val = obj.get_echo;
        end
        function val = get.backscatter(obj)
            val = obj.get_backscatter;
        end
        function val=get.horizontal_position(obj)
            val=obj.horizontal_position_provider.horizontal_position(obj);
        end
        function val=get.vertical_position(obj)
            val=obj.vertical_position_provider.get_vertical_position(obj);
        end
        function val = get.beam_2_instrument_matrix(obj)
            val = obj.instrument_matrix_provider.b2i_matrix(obj);
        end
        function val = get.instrument_2_beam_matrix(obj)
            val = obj.instrument_matrix_provider.i2b_matrix(obj);
        end
        function bad=bad(obj,filter)
            % Mask for profiled data
            %
            %   isbad=bad(obj) returns a mask with ones for data data is marked bad
            %   by the filters.
            %
            %   isbad=bad(obj,filters) allows to specify custom filters, other than
            %   the ones given in the ADCP object.
            %
            %   see also: ADCP, Filter
            if nargin<2
                filter=obj.filters;
            end
            bad=filter(1).bad(obj);
            for co=2:numel(filter)
                bad=bad | filter(co).bad(obj);
            end
        end
        function plot_orientations(obj)
            % Plots orientations of instrument
            %
            %   plot_orientations(obj) plots the pitch, roll and heading.
            %
            % see also: ADCP, plot_velocity, plot_filters, plot_all
            axh=subplot(3,1,1);
            plot(obj.pitch)
            ylabel('Pitch (degr)')
            axh(2)=subplot(3,1,2);
            plot(obj.roll)
            ylabel('Roll (degr)')
            axh(3)=subplot(3,1,3);
            plot(obj.heading)
            ylabel('Head (degr)')
            linkaxes(axh,'x')
        end
        function hfout=plot_velocity(obj,vel)
            % plot the velocity profiles
            %
            %   plot_velocity(obj) plot the velocities in Earth coordinate
            %   system
            %
            %   see also: ADCP
            vel_pos=obj.depth_cell_position;
            vel_pos=mean(vel_pos(:,:,:,3),3,'omitnan');
            if nargin < 2
                vel=obj.velocity(CoordinateSystem.Earth);
            end
            t=obj.time;
            t=seconds(t-t(1));
            hf=gcf;
            axh(1)=subplot(3,1,1);
            pcolor(t,vel_pos,vel(:,:,1));
            hc=colorbar;
            ylabel(hc,'V_x (m/s)')
            shading flat
            axh(2)=subplot(3,1,2);
            pcolor(t,vel_pos,vel(:,:,2));
            ylabel('vertical position (m)')
            hc=colorbar;
            ylabel(hc,'V_y (m/s)')
            shading flat
            axh(3)=subplot(3,1,3);
            pcolor(t,vel_pos,vel(:,:,3));
            hc=colorbar;
            ylabel(hc,'V_z (m/s)')
            shading flat
            xlabel('time (s)')
            linkaxes(axh,'xy')
            if nargout>0
                hfout=hf;
            end
        end
        function hfout=plot_backscatter(obj)
            % plot the backscatter profiles
            %
            %   plot_velocity(obj) plot the backscatter profiles
            %   system
            %
            %   see also: ADCP
            sv_pos=obj.depth_cell_offset;
            sv_pos=mean(sv_pos(:,:,:,3),3,'omitnan');
            sv=obj.backscatter;
            t=obj.time;
            t=seconds(t-t(1));
            nb=size(sv,3);
            hf=gcf;
            axh=nan(nb,1);
            for cb=1:nb
                axh(cb)=subplot(nb,1,cb);
                pcolor(t,sv_pos,sv(:,:,cb));
                clim=mean(sv(:,:,cb),'all','omitnan')+[-2 2]*std(sv(:,:,cb),0,'all','omitnan');
                hc=colorbar;
                ylabel(hc,'Backscatter (dB)')
                shading flat
                set(gca,'clim',clim)
                title(['Beam ',num2str(cb)])
            end
            xlabel('time (s)')
            linkaxes(axh,'xy')
            if nargout>0
                hfout=hf;
            end
        end
        function plot_filters(obj)
            obj.filters.plot(obj);
        end
        function plot_all(obj)
            figure
            obj.plot_orientations;
            figure
            obj.plot_filters;
            figure
            obj.plot_velocity;
            figure
            obj.plot_backscatter;
        end
        function pos=depth_cell_offset(obj,dst)
            % Computes the xyz offset to the profiled depth cells
            %
            %   pos=depth_cell_offset(obj) returns the offset vector (ncells x
            %   nensembles x nbeams x 3) in Earth coordinates
            %
            %   pos = depth_cell_offset(obj,dst) specify the coordinate system
            %   in which the offsets should be returned. dst is a
            %   CoordinateSystem object.
            %
            %   see also: ADCP, depth_cell_position
            if nargin < 2
                dst=CoordinateSystem.Earth;
            end
            tm=-obj.xform(CoordinateSystem.Beam, dst, 'UseTilts', true); % minus since matrix for vel points to adcp
            tm(:,:,:,4)=[];
            pos=tm.*obj.depth_cell_slant_range;
        end
        function pos=depth_cell_position(obj)
            % Computes the xyz positions of the depth cells
            %
            %   pos=depth_cell_position(obj) returns the position vector
            %   (ncells x nensembles x nbeams x 3) in Earth coordinates
            %
            % see also: ADCP, position, depth_cell_offset
            pos=obj.depth_cell_offset + permute([obj.horizontal_position; obj.vertical_position],[3,2,4,1]);
        end
    end
    methods(Abstract, Access=protected)
        val = get_nbeams(obj)
        val = get_nensembles(obj)
        val = get_ncells(obj)
        val = get_coordinate_system(obj)
        val = get_beam_angle(obj)
        val = get_cellsize(obj)
        val = get_time(obj)
        val = get_distmidfirstcell(obj)
        val = get_temperature(obj)
        val = get_salinity(obj)
        val = get_pressure(obj)
        val = get_echo(obj)
        val = get_backscatter(obj)
        val = get_transducer(obj)
    end
    methods(Abstract)
        % Get transformation matrices for coordinate transformations
        %
        %   tm=xform(obj,dst) get the transformation matrices for the
        %   transformation from obj.coordinate_system to the given
        %   destination coordinate system. matrix will be 1 x nensembles x
        %   4 x 4
        %
        %   tm=xform(obj,dst,src) to specify a custom source coordinate
        %   system
        %
        %   tm=xform(...,'ParameterName', parameterValue) allows to
        %   specify the following options:
        %   'UseTilts' - if unset checks the ADCP data in inverse
        %   transformations to know whether to use the tilts. If set to
        %   true will always use the tilts in invers transformations,
        %   if set to false will never use tilts in inverse
        %   transformations.
        %
        %   see also: ADCP
        val = xform(obj,dst,src,varargin)

        % velocity profile data
        %
        %   vel=velocity(obj) returns the profiled velocity in m/s.
        %
        %   vel=velocity(obj,dst) specify destination coordinate system as
        %   CoordinateSystem object.
        %
        %   vel=velocity(obj,dst,filter) specify filter to be used insted
        %   of the ones specified in the current object.
        %
        %   see also: ADCP, CoordinateSystem
        vel = velocity(obj,dst,filter)
    end
end