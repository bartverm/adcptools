classdef Interpolator < handle
    % Base class to define Interpolators
    %
    %   Subclasses should define the interpolate function
    %
    %   obj = Interpolator() default construct Interpolator object
    %
    %   obj = Interpolator(known) construct Interpolator object with known
    %       locations
    %
    % Interpolator properties:
    %   known - the known locations for the interpolation
    %   density_reduction - reduce density of inputs for more efficiency
    %
    % Interpolator properties (read only):
    %   ndims - number of dimensions of the given input points
    %
    % Interpolator methods:
    %   interpolate - interpolate to given query locations (Abstract method)
    %
    %   see also: LoessInterpolator, BathymetryScatteredPoints
    
    properties(SetObservable)
        % Interpolator/known
        %
        %   DxN double array defining the known locations. D is the number of
        %   dimensions + 1 and N the number of points. Defaults to a 3x0 array. The
        %   last row always holds the values to be interpolated
        %
        %   see also: Interpolator
        known double {mustBeFinite}=double.empty(3,0)
        
        % Interpolator/density_reduction
        %
        %   scalar double integer number N that makes the object only use every Nth
        %   input point given. This is usefull when the known locations are very
        %   dense slowing down the computations. Be carefull setting this too high
        %   as it might result in poor interpolations. Default is 1, i.e. use all
        %   known points.
        %
        %   see also: Interpolator
        density_reduction double {mustBeInteger, mustBePositive, mustBeNonzero}=1;
    end
    properties(Dependent, SetAccess=protected, GetAccess=public)
        % Interpolator/n_dims
        %
        %   Returns the number of dimensions of the input points.
        %
        %   see also: Interpolator, known
        n_dims (1,1) double {mustBeFinite,mustBePositive, mustBeInteger}
    end
    methods
        function obj=Interpolator(varargin)
            if nargin > 0
                obj.known=varargin(1);
            end
        end
        function val=get.n_dims(obj)
            val=size(obj.known,1)-1;
        end
        function val=get.known(obj)
            val=obj.known(:,1:obj.density_reduction:end);
        end
        function val=interpolate(obj,query_position) %#ok<STOUT>
            % Interpolates to the given query positions
            %
            %   val=interpolate(obj,query_positions) interpolates to the given
            %   query_positions. The size of query_positions must be DxN with D equal
            %   to the n_dims property and N being the number of query points. The
            %   output has dimensions 1xN.
            %
            %   see also: Interpolator
            validateattributes(query_position,{'numeric'},{'2d','real'},'interpolate','query_position',2)
            assert(size(query_position,1)==obj.n_dims,['query_position must have ',num2str(obj.n_dims),' dimensions'])
        end
    end
end