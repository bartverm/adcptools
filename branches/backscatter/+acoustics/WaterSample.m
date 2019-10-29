classdef WaterSample < handle
% acoustics.WaterSample water sample for sediment concentration estimates
%
%   acoustics.WaterSample properties:
%   distance - distance from transducer to the water sample
%   sediment_density - density of sediment, or particles in water
%   time - time the sample was collected
%   distribution - grain size distribution of the sample
%   volume_concentration - volume concentration of sediment or particles
%   mass_concentration - mass concentration of sample
%
%   acoustics.WaterSample methods:
%   backscatter_strength - estimate backscatter strength of sample
%
%   see also: acoustics.GrainSizeDistribution
    properties
        % acoustics.WaterSample/distance
        %
        % distance to the water sample from the transducer. Default is 0
        %
        %  see also:acoustics.WaterSample
        distance(1,1) double {mustBeFinite, mustBeNonnegative} = 0
        
        % acoustics.WaterSample/mass_concentration
        %
        % mass concentration of the water sample. Value must be a scalar,
        % non-negative double. Default is 0;
        %
        % see also: acoustics.WaterSample
        mass_concentration(1,1) double {mustBeFinite, mustBeNonnegative} = 0

        % acoustics.WaterSample/volume_concentration
        %
        % volume concentration of the water sample. Value must be a scalar,
        % non-negative double. Default is 0;
        %
        % see also: acoustics.WaterSample
        volume_concentration(1,1) double {mustBeFinite, mustBeNonnegative} = 0

        % acoustics.WaterSample/time
        %
        % time the sample was taken. This is a scalar datetime object.
        %
        % see also: acoustics.WaterSample, datetime
        time(1,1) datetime = datetime()
        
        % acoustics.WaterSample/distribution
        %
        % grain size distribution of the particles in the sample. This is a
        % scalar acoustics.GrainSizeDistribution object. Default is an
        % empty object.
        %
        % see also: acoustics.WaterSample, acoustics.GrainSizeDistribution
        distribution(:,1) acoustics.GrainSizeDistribution = acoustics.GrainSizeDistribution.empty
    end
    properties(Dependent, SetAccess=private)
        % acoustics.WaterSample/sediment_density
        %
        % density of the (sediment) particles in the sample. This is a
        % read-only property computed from mass and sediment concentrations
        %
        % see also: acoustics.WaterSample, mass_concentration,
        % volume_concentration
        sediment_density(1,1) double {mustBeFinite, mustBeNonnegative}
    end
    methods % set and get methods
        function set.distribution(obj,val)
            assert(isscalar(val),'Property distribution must be scalar')
            obj.distribution=val;
        end
        function val=get.sediment_density(obj)
            val=obj.mass_concentration./obj.volume_concentration+1000;
        end
    end
    methods % ordinary methods
        function sv=backscatter_strength(obj,wavenumber)
        % Compute the backscatter strength of the water sample
        %
        %   sv=backscatter_strength(obj,wavenumber) computes the volume
        %   backsactter strength of the water sample for the given
        %   wavenumber. Wavenumber can also be a PistonTransducer object.
        %   This function requires the distribution property to be
        %   non-empty.
        %
        %   see also: acoustics.WaterSample, distribution
            obj.distrib_not_empty('compute backscatter strength');
            sv=10*log10(obj.distribution.ks_squared(wavenumber)*obj.mass_concentration);
        end
            
    end
    methods(Access=protected)
        function distrib_not_empty(obj,what)
        % generates an error if the distribution property is empty
            if nargin < 2
                what='complete operation';
            end
            assert(~isempty(obj.distribution),'WaterSample:empty_distribution',['''distribution'' property cannot be empty to ', what])
        end
    end
end