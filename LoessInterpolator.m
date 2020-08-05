classdef LoessInterpolator < Interpolator
    properties
        span(1,1) double {mustBePositive, mustBeFinite} = 0.1;
        robust_iterations(1,1) double {mustBeNonnegative, mustBeInteger}=0;
        n_threads(1,1) double {mustBePositive, mustBeFinite, mustBeInteger}=maxNumCompThreads;
        order(1,1) double {mustBePositive, mustBeFinite,mustBeInteger,mustBeLessThan(order,3)}=1;
    end
    methods
        function val=interpolate(obj,query_position)
            validateattributes(query_position,{'numeric'},{'2d','real'},'interpolate','query_position',2)
            assert(size(query_position,1)==obj.n_dims,['query_position must have ',num2str(obj.n_dims),' dimensions'])
            
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