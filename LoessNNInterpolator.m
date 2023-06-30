classdef LoessNNInterpolator < Interpolator
% Loess, Natural neighbour interpolator class
%
%   Performs smoothing through robust and local weighted regressions. This 
%   class uses the C++ loess function (see 
%   https://github.com/bartverm/loess). The class smooths points and
%   filters out outliers with loess, uses natural neighbour interpolation
%   and nearest neighbour extrapolation. It stores the internally the
%   smoothed points in a scatteredInterpolant, which can be reused, until
%   the 'known' property is modified, invalidating the internal
%   interpolator.
%
%   LoessNNInterpolator properties:
%   known - known positions and values for the interpolation
%   span - defines the fraction of points to perform the local regressions
%   robust_iterations - number of robust iterations to filter out outliers
%   n_threads - number of computations threads to use
%   order - order of the regression
%
%   LoessNNInterpolator properties (read only):
%   n_dims - number of dimensions of known points
%
%   LoessNNInterpolator methods:
%   interpolate - performs the interpolation
%
%   see also: Interpolator, GridddataInterpolator, LoessInterpolator
    properties(SetObservable=true)

% LoessInterpolator/span
%
%   scalar double (0 > span <= 1) defining the fraction of points to be
%   used for the local regression. The larger this number the smoother the
%   result, but the slower the computation. Default value is 0.01, i.e. 1%
%   of the input points are included in the local regression. If the number
%   is >1, the span is interpreted as number of points.
%
% see also: LoessInterpolator
        span(1,1) double {mustBePositive, mustBeFinite} = 0.01;

% LoessInterpolator/robust_iterations
%
%   scalar double indicating the number of robust iterations that should be
%   performed to exclude outliers. Each robust iteration will remove
%   outliers based on the residuals. Points exceeding six times the median
%   residual will be excluded from the regression. Two robust iteration
%   will remove most outliers. Every robust iterations will slow the
%   interpolation. Default is 0, so no outlier removal
%
%   see also: LoessInterpolator
        robust_iterations(1,1) double {mustBeNonnegative, mustBeInteger}=0;
        
% LoessInterpolator/n_threads
%
%   Number of computational threads to be used for te interpolation. This
%   defaults to maxNumCompThreads. This can be changed up to the maximum
%   number of computational threads available on the system. Note that
%   initializing a parallel pool on matlab sometimes slows down the
%   multithreading in the loess function.
%
% see also: LoessInterpolator
        n_threads(1,1) double {mustBePositive, mustBeFinite, mustBeInteger}=maxNumCompThreads;
        
% LoessInterpolator/order
%
%   This can be either 1 or 2 and defines the order of the local
%   regression.
%
% see also: LoessInterpolator
        order(1,1) double {mustBePositive, mustBeFinite,mustBeInteger,mustBeLessThan(order,3)}=1;
    end
    properties(Access=protected)
        int(1,1) scatteredInterpolant
        update_int(1,1) logical
    end
    methods
        function obj=LoessNNInterpolator(varargin) % constructor
            assert(helpers.load_loess_package,...
                'Failed to load the loess submodule')
            obj@Interpolator(varargin{:}) % call parent constructor
            obj.update_int=true; % flag the interpolator to be updated
            addlistener(obj,'known','PostSet',@obj.reset_interpolant); % create property listener to trigger new smoothing
            addlistener(obj,'span','PostSet',@obj.reset_interpolant); % create property listener to trigger new smoothing
            addlistener(obj,'robust_iterations','PostSet',@obj.reset_interpolant); % create property listener to trigger new smoothing
            addlistener(obj,'order','PostSet',@obj.reset_interpolant); % create property listener to trigger new smoothing
        end
        function reset_interpolant(obj,varargin)
            obj.update_int=true; % makes sure smoother is called upon next interpolation
        end
        function make_interpolant(obj) % Smooth points and create the interpolant

            val=loess(obj.known(1:end-1,:)',... % call the loess function
                obj.known(end,:)',...
                obj.known(1:end-1,:)',...
                obj.span,...
                obj.robust_iterations,...
                obj.order,...
                obj.n_threads);           
            obj.int=scatteredInterpolant(obj.known(1:end-1,:)', val,'natural','nearest'); % create the interpolant
            obj.update_int=false; % make sure smoother is not called again when interpolating
        end
        function val=interpolate(obj,query_position)
            if obj.update_int 
                obj.make_interpolant()
            end
            interpolate@Interpolator(obj,query_position); % call parent interpolate function
            val=obj.int(query_position'); % perform interpolation
        end
    end
end