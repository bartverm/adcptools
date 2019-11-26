classdef Water < handle
    % acoustics.Water defines typical water properties
    %
    %   acoustics.Water properties:
    %       temperature - Temparature of water in Celsius
    %       salinity - Salinity of water in ppm
    %
    %   acoustics.Water dependent properties (read only)
    %       density - Density of water in kg/m^3
    %       dyn_viscosity - Dynamic viscosity of water kg/m/s
    %       kin_viscosity - Kinematic viscosity of water m^2/s
    %
    %   acoustics.Water methods:
    %       speedsound - Computes speed of sound in water
    %   
    %   acoustics.Water static methods:
    %       speedsound_static - Computes speed of sound in water
    %       density_static - Compute density of water
    %       dyn_viscosity_static - Compute dynamic viscosity of water
    %       kin_viscosity_static - Compute kinematic viscosity of water
    %
    % see also: acoustics
    
    properties(SetObservable, AbortSet)
        temperature double {mustBeFinite}=20 % Temperature of water in Celsius, default is 20
        salinity double {mustBeNonnegative}=200  % Salinity of water in ppm (parts per million), default is 200
        pH double {mustBeFinite}=7 % pH of water
    end
    properties(Dependent)
        % Density of seawater
        %
        % Computation according to equation 8 in <a href="matlab:web('https://dx.doi.org/10.5004/dwt.2010.1079')">Sharqawy et al., (2012)</a>.
        density
        
        % Dynamic viscosity of seawater
        %
        % Computation accaccording to equation 22 and 23 in <a href="matlab:web('https://dx.doi.org/10.5004/dwt.2010.1079')">Sharqawy et al., (2012)</a>.
        dyn_viscosity
        
        % Kinematic viscosity of seawater
        %
        % Computed from dynamic viscosity and density
        kin_viscosity
    end
    methods
       % Property get methods for public properties
        function val=get.density(obj)
            try
                val=acoustics.Water.density_static(obj.temperature, obj.salinity);
            catch err
                acoustics.Water.handle_err(err);
            end
        end
        function val = get.dyn_viscosity(obj)
            try
                val=acoustics.Water.dyn_viscosity_static(obj.temperature, obj.salinity);
            catch err
                acoustics.Water.handle_err(err)
            end
        end
        function val = get.kin_viscosity(obj)
            try
                val=acoustics.Water.kin_viscosity_static(obj.temperature, obj.salinity);
            catch err
                acoustics.Water.handle_err(err)
            end
        end
        
        % Regular methods
        function value=speedsound(obj,depth)
        % Computes speed of sound in water
        %
        %   C=obj.speedsound(DEPTH) computes speed of sound C (m/s) in
        %   water at the given DEPTH (m). Uses Medwin and Clay, pg 85, Eq 3.3.3
        %
        %   see also: acoustics.Water.speedsound_static
            value=acoustics.Water.speedsound_static(obj.temperature,obj.salinity,depth);
        end
    end
    methods (Static,Access=private)
        function handle_err(err)
            if strcmp(err.identifier,'elementwise_output_size:incompatible')
                error('property dimensions are incompatible for elementwise computations');
            end
        end
    end
    methods (Static)
        function c=speedsound_static(temperature, salinity, depth)
        % Computes speed of sound in water
        %
        %   C=acoustics.Water.speedsound(T,S,D) computes speed of sound C (m/s) in
        %   water at the given temperature T (celsius), salinity S (ppm) depth D (m).
        %   Uses Medwin and Clay, pg 85, Eq 3.3.3
        %
        %   see also: acoustics.Water.speedsound
           helpers.elementwise_output_size(temperature, salinity, depth);
           c = 1448.96 + 4.591     * temperature-...
                         5.304e-2  * temperature.^2+...
                         2.374e-4  * temperature.^3+...
                         1.340     * (salinity/1000-35)+...
                         1.630e-2  * depth+1.675e-7*depth.^2-...
                         1.025e-2  * temperature.*(salinity/1000-35)-...
                         7.139e-13 * temperature.*depth.^3;                %Speed of sound, Medwin and Clay, pg 85, Eq 3.3.3
          % TODO: make values NaN outside of calibration range!
        end
        function rho=density_static(temperature, salinity)
        % Computes the density of saltwater
        %
        %   RHO = acoustics.Water.density_static(TEMPERATURE, SALINITY)
        %   computes the density (kg/m^3) of saltwater given TEMPERATURE (celsius) and
        %   SALINITY (ppm). Array dimensions must be compatible. 
        %   Computations based on equation 8 in <a href="matlab:web('https://dx.doi.org/10.5004/dwt.2010.1079')">Sharqawy et al., (2012)</a>.
        %   Valid range 0 < t < 180, 0 < S < 1.6e6. Outside of this range
        %   density will be NaN
        %
        %   see also: acoustics.Water.density
            helpers.elementwise_output_size(temperature, salinity)
            t=temperature;
            S=salinity/1e6;
            a1 = 9.999e2;
            a2 = 2.034e-2;
            a3 = -6.162e-3;
            a4 = 2.261e-5;
            a5 = -4.657e-8;
            b1 = 8.020e2;
            b2 = -2.001;
            b3 = 1.677e-2;
            b4 = -3.060e-5;
            b5 = -1.613e-5;
            rho =   (a1 + a2 * t + a3 * t .^ 2 + a4 * t .^ 3 + a5 * t .^ 5) +...
                S .* (b1 + b2 * t + b3 * t .^ 2 + b4 * t .^ 3 + b5 * S .^ 2 .* t .^ 2);
            rho(t < 0 | t> 180 | S < 0 | S > 0.16) = nan; % out of calibration range
        end
        function mu=dyn_viscosity_static(temperature, salinity)
        % Computes the dynamic viscosity of salt water
        %
        %   RHO = acoustics.Water.dyn_viscosity_static(TEMPERATURE, SALINITY)
        %   computes the dynamic viscosity (kg/m/s) of saltwater given TEMPERATURE (celsius) and
        %   SALINITY (ppm). Array dimensions must be compatible. 
        %   Computations based on equation 22 and 23 in <a href="matlab:web('https://dx.doi.org/10.5004/dwt.2010.1079')">Sharqawy et al., (2012)</a>.
        %   Valid range 0 < t < 180, 0 < S < 1.5e6. Outside of this range
        %   result will be NaN
        %
        %   see also: acoustics.Water.dyn_viscosity
            helpers.elementwise_output_size(temperature, salinity)
            t=temperature;
            S=salinity/1e6;
            mu_w = 4.2844e-5 + (0.157 * ( t + 64.993 ) .^ 2 - 91.296 ) .^ (-1);
            A = 1.541 + 1.998e-2 * t - 9.52e-5 * t .^ 2;
            B = 7.974 - 7.561e-2 * t + 4.724e-4 * t .^ 2;
            mu = mu_w .* (1 + A .* S + B .* S .^ 2);
            mu(t < 0 | t > 180 | S < 0 | S > 0.15)=nan;
        end
        function nu=kin_viscosity_static(temperature, salinity)
        % Computes the kinematic viscosity of salt water
        %
        %   RHO = acoustics.Water.kin_viscosity_static(TEMPERATURE, SALINITY)
        %   computes the kinematric viscosity (m^2/s) of saltwater given TEMPERATURE (celsius) and
        %   SALINITY (ppm). Array dimensions must be compatible. 
        %   Computations based on equation 8, 22 and 23 in <a href="matlab:web('https://dx.doi.org/10.5004/dwt.2010.1079')">Sharqawy et al., (2012)</a>.
        %   Valid range 0 < t < 180, 0 < S < 1.5e6. Outside of this range
        %   result will be NaN
        %
        %   see also: acoustics.Water.kin_viscosity
            nu = acoustics.Water.dyn_viscosity_static(temperature,salinity)./...
                 acoustics.Water.density_static(temperature, salinity);
        end
    end
end