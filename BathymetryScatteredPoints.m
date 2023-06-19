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
        known (3,:) double {mustBeFinite} = double.empty(3, 0);

        % BathymetryScatteredPoints/interpolator performs interpolation
        %
        %   Scalar Interpolator object performing the interpolation. Default is a
        %   LoessInterpolator object.
        %
        %   see also: BathymetryScatteredPoints, LoessInterpolator,
        %   LoessNNinterpolator
        interpolator (1,1) Interpolator = LoessNNInterpolator;
    end
    methods
        function obj=BathymetryScatteredPoints(varargin)
            obj = obj@Bathymetry(varargin{:});
            all_arg_in = varargin;
            varargin = obj(1).unprocessed_construction_inputs;
            nargin = numel(varargin);
            isprocessed = true(size(varargin));
            siz = size(obj);
            siz = num2cell(siz);
            int(siz{:}) = LoessNNInterpolator;
            % make sure interpolators are different for array
            obj.assign_property('interpolator', int);
            addlistener(obj, 'known',...
                'PostSet',@obj.set_interpolator_known);
            addlistener(obj,'interpolator',...
                'PostSet',@obj.set_interpolator_known);
            construct_from_vmadcp = false;
            construct_water_level = false;
            has_filt = false;
            has_vmadcp = false;
            for ca = 1 : nargin
                cur_arg = varargin{ca};
                if isa(cur_arg, 'VMADCP')
                    construct_from_vmadcp = true;
                    construct_water_level = true;
                    vadcp = cur_arg;
                    has_vmadcp = true;
                    continue
                elseif isa(cur_arg, 'EnsembleFilter')
                    filter = cur_arg;
                    has_filt = true;
                    continue
                elseif isa(cur_arg, 'Interpolator')
                    var_name = 'interpolator';
                elseif isa(cur_arg, 'WaterLevel')
                    construct_water_level = false;
                    continue
                elseif isa(cur_arg, 'double')
                    var_name = 'known';
                else
                    isprocessed(ca) = false;
                    continue
                end
                obj.assign_var(var_name, cur_arg)
            end
            varargin(isprocessed) = [];
            obj(1).unprocessed_construction_inputs = varargin;

            if construct_from_vmadcp && has_vmadcp
                if ~has_filt
                    f_exp = find(strcmp('NoExpand',all_arg_in));
                    f_adcp = find(cellfun(@(x) isa(x, 'VMADCP'), all_arg_in));
                    args_idx = sort([f_exp f_adcp]);
                    filter = EnsembleFilter(all_arg_in{args_idx});
                end
                obj.known_from_vmadcp(vadcp, filter)
            end
            if construct_water_level
                wls = {vadcp.water_level_object};
                if isscalar(obj) && ~isscalar(wls)
                    if ~isequal(wls{:})
                        warning(['Using water level from', ...
                            ' first ADCP object']);
                    end
                    wls = wls(1);
                end
                obj.assign_property('water_level', [wls{:}]);
            end
        end
        function z = get_bed_elev(obj, x, y)
            validateattributes(x, {'double'}, {})
            validateattributes(y, {'double'}, {})
            assert(isequal(size(x), size(y)),...
                'Size of x and y should match')
            sizin = size(x);
            z = reshape(obj.interpolator.interpolate([...
                reshape(x ,1, []);...
                reshape(y,1,[])]),sizin);
        end
        function known_from_vmadcp(obj, varargin)
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
            if ~isscalar(obj)
                obj.run_method('known_from_vmadcp', varargin{:});
                return
            end
            vmadcp = [];
            construct_filter = true;
            for ci = 1:numel(varargin)
                carg = varargin{ci};
                if isa(carg, 'VMADCP')
                    vmadcp = carg;
                elseif isa(carg,'EnsembleFilter')
                    filter = carg;
                    construct_filter=false;
                end
            end
            assert(~isempty(vmadcp),'A VMADCP object is required as input')
            if construct_filter
                filter=EnsembleFilter(vmadcp);
            end
            tpos = vmadcp.cat_property('bed_position');
            xpos = tpos(:, :, :, 1);
            ypos = tpos(:, :, :, 2);
            zpos = tpos(:, :, :, 3);
            isfin = all(isfinite(tpos), 4);
            var = cell(size(filter));
            for cc = 1 : numel(var)
                filt = isfin & ~filter(cc).all_cells_bad(vmadcp);
                var{cc} = [xpos(filt)'; ypos(filt)'; zpos(filt)'];
            end
            obj.assign_property('known', var);
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
            z_interp = obj.get_bed_elev(obj.known(1, :), obj.known(2, :));
            res = z_interp - obj.known(3, :);
            scatter(obj.known(1, :), obj.known(2, :), 5, res, 'filled')
            set(gca, 'clim', mean(res, 'omitnan') + ...
                std(res, 'omitnan') * 2 * [-1 1])
            colorbar
            axis equal
        end
        function varargout = plot(obj)
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
            hp=plot3(obj.known(1,:),obj.known(2,:),obj.known(3,:),...
                'k.','MarkerSize',.1);
            hold on
            pbaspect([5 5 1])
            da = daspect;
            hrat = max(da(1:2))/da(3);
            daspect([hrat hrat 1])
            as=alphaShape(obj.known(1,:)',obj.known(2,:)',1);
            as.Alpha=as.criticalAlpha('one-region')*6;
            tri = alphaTriangulation(as);
            z_interp=obj.get_bed_elev(obj.known(1,:),obj.known(2,:));
            ht=trimesh(tri,obj.known(1,:), obj.known(2,:), z_interp ,...
                'FaceColor','interp','EdgeColor','none');
            set(gca,'NextPlot',hold_stat)
            set(gca,'Clipping','off')
            if nargout > 0
                varargout{1} = hp;
            end
            if nargout > 1
                varargout{2} = ht;
            end
        end % function
    end % methods
    methods (Static, Access=protected)
        function set_interpolator_known(varargin)
            srcObj = varargin{2}.AffectedObject;
            srcObj.interpolator.known = srcObj.known;
        end % function
    end % protected methods
end % classdef