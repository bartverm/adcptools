classdef LoessInterpolator < Interpolator
% Loess interpolator class
%
%   Performs interpolations and smoothing through robust and local 
%   weighted regressions. This class is a wrapper of the C++ loess function
%   (see https://github.com/bartverm/loess). The class will also
%   extrapolate when points are outside the convex hull of the known
%   inputs.
%
%   LoessInterpolator properties:
%   known - known positions and values for the interpolation
%   span - defines the fraction of points to perform the local regressions
%   robust_iterations - number of robust iterations to filter out outliers
%   n_threads - number of computations threads to use
%   order - order of the regression
%
%   LoessInterpolator properties (read only):
%   n_dims - number of dimensions of known points
%
%   LoessInterpolator methods:
%   interpolate - performs the interpolation
%
%   see also: Interpolator, GridddataInterpolator
    properties

% LoessInterpolator/span
%
%   scalar double (0 > span <= 1) defining the fraction of points to be
%   used for the local regression. The larger this number the smoother the
%   result, but the slower the computation. Default value is 0.01, i.e. 1%
%   of the input points are included in the local regression
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
%   interpolation.
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
    methods
        function val=interpolate(obj,query_position)
            interpolate@Interpolator(obj,query_position);
            
            % call to loess function
            fl=exist('loess','file');
            if ~(fl==2 || fl==3)
                error('could not find loess function, please add it to the path')
            end
            try
                val=loess(obj.known(1:end-1,:)',...
                    obj.known(end,:)',...
                    query_position',...
                    obj.span,...
                    obj.robust_iterations,...
                    obj.order,...
                    obj.n_threads);
            catch err
                throw(err)
            end
        end
    end
end