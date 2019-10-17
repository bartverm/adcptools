classdef WaterSample < handle
    properties
        distance(1,1) double {mustBeFinite, mustBeNonnegative} = 0
        mass_concentration(1,1) double {mustBeFinite, mustBeNonnegative} = 0
        sediment_density(1,1) double {mustBeFinite, mustBeNonnegative} = 2650
        time(1,1) datetime = datetime()
        distribution(:,1) acoustics.GrainSizeDistribution = acoustics.GrainSizeDistribution.empty
    end
    properties(Dependent)
        volume_concentration(1,1) double {mustBeFinite, mustBeNonnegative}
    end
    methods
        function set.distribution(obj,val)
            assert(isscalar(val),'Property distribution must be scalar')
            obj.distribution=val;
        end
        function val=get.volume_concentration(obj)
            val=obj.mass_concentration/obj.sediment_density;
        end
        function set.volume_concentration(obj,val)
            obj.mass_concentration=val*obj.sediment_density;
        end
        function n=n_particles(obj)
            obj.distrib_not_empty('compute the number of particles');
            n=obj.distribution.n_particles(obj.volume_concentration);
        end
        function sv=backscatter_strength(obj,wavenumber)
            obj.distrib_not_empty('compute backscatter strength');
            sv=10*log10(obj.distribution.bs_xsection(wavenumber)*obj.n_particles);
        end
            
    end
    methods(Access=protected)
        function distrib_not_empty(obj,what)
            if nargin < 2
                what='complete operation';
            end
            assert(~isempty(obj.distribution),'WaterSample:empty_distribution',['''distribution'' property cannot be empty to ', what])
        end
    end
end