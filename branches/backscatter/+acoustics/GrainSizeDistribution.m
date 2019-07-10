classdef GrainSizeDistribution < handle
    properties
        grainsize(:,1) double {mustBeNonnegative, mustBeNonNan, mustBeNonempty} = 0
        fraction(:,1) double {mustBeFinite, mustBeNonempty, mustBeGreaterThanOrEqual(fraction,0), mustBeLessThanOrEqual(fraction, 1)}=0
    end
    properties(Dependent)
        cumulative_fraction(:,1) double {mustBeFinite, mustBeNonempty, mustBeGreaterThanOrEqual(cumulative_fraction,0) mustBeLessThanOrEqual(cumulative_fraction, 1)}
    end
    methods
        function set.cumulative_fraction(obj,val)
            obj.fraction=diff(val);
        end
        function val=get.cumulative_fraction(obj)
            val=cumsum(obj.fraction);
        end
        function val=D(N)
            obj.check_grainsize_fraction();
            
        end
        function check_grainsize_fraction(obj)
            assert(isequal(size(obj.grainsize), size(obj.fraction)),'fraction and grainsize properties should have the same size')
            assert(sum(obj.fraction)-1 < eps,'sum of all fractions should equal 1')
            assert(issorted(obj.grainsize),'grainsize must be in ascending order')
        end
        function plot(obj)
            plot(obj.grainsize,obj.fraction)
            set(gca,'xscale','logarithmic')
        end
    end
end