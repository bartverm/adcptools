classdef BathymetryScatteredPoints < Bathymetry
% creates bathymetry from scattered input points
%
%   This class is a concrete implementation of the Bathymetry class
%
%   the bathymetry is obtained interpolating the bathymetry from scattered
%   input points defined in the property 'known' using the Interpolator
%   object defined in the property 'interpolator'. 
%
%   obj=BathymetryScatteredPoints(...) uses objects passed on construction
%   to assign the object properties according to the class of the arguments
%   passed:
%   - VMADCP objects are used to assign ADCP bed detections as input points
%   - EnsembleFilter is used to only include ADCP points defined with the
%       filter
%   - Interpolator object is assigned to the 'interpolator' property
%   - WaterLevel objects are assigned to the 'water_level' property
%
%   If a VMADCP object is passed upon construction and no WaterLevel object
%   is given, the water_level property is set from the VMADCP object
%
%   BathymetryScatteredPoints properties:
%   known - known scatterd input points (3xN)
%   interpolator - Interpolator object performing the interpolation
%
%   BathymetryScatteredPoints methods:
%   plot - plots the fitted bathymetry and the scattered input points
%   plot_residuals - plots the residuals of the interpolation at the inputs
%   get_bed_elev - returns bed elevation at given query points
%   get_depth - returns depth at given query points at given time
%   known_from_vmadcp - sets the known points from VMADCP object
%
%   see also: Bathymetry, Interpolator, VMADCP, EnsembleFilter, WaterLevel
    properties (SetObservable)
% BathymetryScatteredPoints/known known points for the interpolation
%
%   Known points of the interpolation given as a 3xN double array holding
%   the x,y,z coordinates of the bed for N points. Only finite points are
%   allowed
%
%   see also: BathymetryScatteredPoints
        known (3,:) double {mustBeFinite} = double.empty(3,0);

% BathymetryScatteredPoints/interpolator performs interpolation
%
%   Scalar Interpolator object performing the interpolation. Default is a
%   LoessInterpolator object.
%
%   see also: BathymetryScatteredPoints, LoessInterpolator,
%   LoessNNinterpolator
        interpolator (1,1) Interpolator;
    end
    methods
        function obj=BathymetryScatteredPoints(varargin)
            obj=obj@Bathymetry(varargin{:});
            obj.interpolator=LoessNNInterpolator;
            addlistener(obj,'known','PostSet',@obj.set_interpolator_known);
            addlistener(obj,'interpolator','PostSet',@obj.set_interpolator_known);
            filter=EnsembleFilter;
            construct_from_vmadcp=false;
            construct_water_level=false;
            for ca=1:nargin
                cur_arg=varargin{ca};
                if isa(cur_arg,'VMADCP')
                    construct_from_vmadcp=true;
                    construct_water_level=true;
                    vadcp=cur_arg;
                elseif isa(cur_arg,'EnsembleFilter')
                    filter=[filter, cur_arg]; %#ok<AGROW>
                elseif isa(cur_arg,'Interpolator')
                    obj.interpolator=cur_arg;
                elseif isa(cur_arg,'WaterLevel')
                    construct_water_level=false;
                elseif isa(cur_arg,'double')
                    obj.known = cur_arg;
                else
                    warning('Bathymetry:unhadled_input',['Unhandled input of type: ', class(cur_arg)])
                end
            end
            if construct_from_vmadcp
                obj.known_from_vmadcp(vadcp,filter)
            end
            if construct_water_level
                obj.water_level=vadcp.water_level_object;
            end
        end
        function z=get_bed_elev(obj,x,y)
            validateattributes(x,{'double'},{})
            validateattributes(y,{'double'},{})
            assert(isequal(size(x),size(y)),'Size of x and y should match')
            sizin=size(x);
            z=reshape(obj.interpolator.interpolate([reshape(x,1,[]); reshape(y,1,[])]),sizin);
        end
        function known_from_vmadcp(obj,vmadcp,filter)
% BathymetryScatteredPoints/known_from_adcp
%
%   obj.known_from_vmadcp(vmadcp) sets the known scattered points
%   for the interpolation from the bed detections of the vessem mounted
%   ADCP.
%
%   obj.known_from_vmadcp(vmadcp, filter) allows to set a filter to exclude
%   part of the ensembles.
%
%   see also: BathymetryScatteredPoints, VMADCP, VMADCP/bed_position
            validateattributes(vmadcp,{'VMADCP'},{'scalar'},'pos_from_vmadcp','vmadcp',2)
            tpos=vmadcp.bed_position;
            tpos(:,filter.all_cells_bad(vmadcp),:,:)=nan;
            xpos=tpos(:,:,:,1);
            ypos=tpos(:,:,:,2);
            zpos=tpos(:,:,:,3);
            isfin=all(isfinite(tpos),4);
            obj.known=[xpos(isfin)';ypos(isfin)';zpos(isfin)'];
        end
        function plot_residuals(obj)
% Plots the residuals at the input points
%
%   obj.plot_residuals, plots the residuals of the smoothing at the input
%   points. It is important that the residuals are more or less randomly
%   distributed without obvious spatial structure. If the residuals feature
%   spatial structure (e.g. zones whith residuals always positive) consider
%   reducing the level of smoothing. Colors are scaled around the mean
%   residuals +/- 2 standard deviations.
%
%   see also: BathymetryScatteredPoints
            z_interp=obj.get_bed_elev(obj.known(1,:),obj.known(2,:));
            res=z_interp-obj.known(3,:);
            scatter(obj.known(1,:),obj.known(2,:),5,res,'filled')
            set(gca,'clim',mean(res,'omitnan') + std(res,'omitnan')*2*[-1 1])
            colorbar
            axis equal
        end
        function varargout=plot(obj)
% plots the bathymetry
%
%   obj.plot() creates a 3D plot of the interpolated bathymetry surface.
%   The surface is a triangulated surface between the input location. The
%   surface does usually not cross the input points due to the smoothing in
%   the interpolator. The surface extends in the alphaShape of the input
%   points, with an alphaShape radius of 6 times the radius giving one
%   output region. 
%
%   It also plots the input points as black dots.
%
%   [hp, ht]=obj.plot returns a handle to the plotted input points, and to
%   the triangulated surface representing the bathymetry.
%
%   obj.plot can also plot an array of objects which will include several
%   bathymetries in the plot. This feature is handled by the Bathymetry
%   superclass.
%
%   see also: BathymetryScatteredPoints, plot_residuals, Bathymetry
           if ~isscalar(obj)
               plot@Bathymetry(obj);
               return
           end
           hold_stat=get(gca,'NextPlot');
           hp=plot3(obj.known(1,:),obj.known(2,:),obj.known(3,:),'k.','MarkerSize',.1);
           hold on
           set(gca,'dataaspectratio',[5 5 1])
           as=alphaShape(obj.known(1,:)',obj.known(2,:)',1);
           as.Alpha=as.criticalAlpha('one-region')*6;
           tri = alphaTriangulation(as);
           z_interp=obj.get_bed_elev(obj.known(1,:),obj.known(2,:));
           ht=trimesh(tri,obj.known(1,:), obj.known(2,:), z_interp ,'FaceColor','interp','EdgeColor','none');
           set(gca,'NextPlot',hold_stat)
           if nargout > 0
               varargout{1}=hp;
           end
           if nargout > 1
               varargout{2}=ht;
           end
        end % function
    end % methods
    methods (Access=protected)
        function set_interpolator_known(obj,varargin)
            obj.interpolator.known=obj.known;
        end % function
    end % protected methods
end % classdef