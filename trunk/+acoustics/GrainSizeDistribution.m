classdef GrainSizeDistribution < handle
    % acoustics.GrainSizeDistribution computes acoustic sediment properties
    %
    %   obj=GrainSizeDistribution() Construct a default object
    %
    %   obj=GrainSizeDistribution(diameter) specify the diameter property
    % 
    %   obj=GrainSizeDistribution(diameter, fraction) specify both the 
    %       diameter and the fraction properties
    %
    %   acoustics.GrainSizeDistribution properties:
    %       diameter - diameters of class boundaries
    %       fraction - fractions for each class
    %       radius - radius of class boundaries
    %       cumulative_fraction - cumulative fractions for each class
    %
    %   acoustics.GrainSizeDistribution methods:
    %       D - get the N-th percentile of the distribution
    %       gsd_int - integrate a function over the distribution
    %       moment - compute moments of the distribution
    %       mean_diam - compute the mean diameter
    %       mean_radius - compute the mean radius
    %       form_function - compute the form function values
    %       form_function_static - static version of above function
    %       norm_scatter_xsection - normalized scattering cross-section
    %       norm_scatter_xsection_static - static version of above function
    %       bs_xsection - mean backscattering cross section
    %       ks_squared - ks_squared for volume backscattering computation
    %       specific_atten - compute mean specific attenuation
    %       check_diameter_fraction - check properties
    %       plot - plot the GSD
    %       plot_hist - plot the GSD as a histogram
    %       plot_cumulative - plot the cumulative GSD
    %       plot_form_function - plot the form function
    %       plot_all - plot the GSD and the form function
    %
    %   see also: acoustics.WaterSample, acoustics.Water
    
    
    properties
        % acoustics.GrainSizeDistribution/diameter property
        %
        %   Vector with grain diameters of the class boundaries. This
        %   vector should be one element longer than the fraction and
        %   cumulative fraction property. This property is linked with the
        %   radius property. Any changes to the radius property will also
        %   affect the diameters. Values should not be NaN, or empty, they
        %   must be non-negative and sorted. Default is [0; inf]
        %
        %   see also: acoustics.GrainSizeDistribution, radius
        diameter(:,1) double {mustBeNonnegative, mustBeNonNan, mustBeNonempty} = [0; inf]
        
        % acoustics.GrainSizeDistribution/fraction property
        %
        %   Vector with fractions of grains that are within a certain
        %   class. The vector must have one element less than the diameter 
        %   property. This property is linked to the cumulative fraction
        %   property. Altering one will also alter the other. Values must
        %   be finite, non-empty, larger or equal to 0 and smaller or equal
        %   to 1. The sum of all fractions should be equal to one. Default
        %   is 1;
        %
        %   see also: acoustics.GrainSizeDistribution, cumulative_fraction
        fraction(:,1) double {mustBeFinite, mustBeNonempty, mustBeGreaterThanOrEqual(fraction,0), mustBeLessThanOrEqual(fraction, 1)}=1
    end
    properties(Dependent)
        % acoustics.GrainSizeDistribution/radius property
        %
        %   Vector with grain radiuses of the class boundaries. This
        %   vector should be one element longer than the fraction and
        %   cumulative fraction property. This property is linked with the
        %   diameter property. Altering one will also affect the other. 
        %   Values should not be NaN, or empty, they must be non-negative
        %   and sorted. Default is [0; Inf]
        %
        %   see also: acoustics.GrainSizeDistribution, diameter
        radius(:,1) double {mustBeNonnegative, mustBeNonNan, mustBeNonempty}
        
        % acoustics.GrainSizeDistribution/cumulative_fraction property
        %
        %   Vector with cumulative fractions of grains that are within a 
        %   certain class. The vector must have one element less than the 
        %   diameter or radius property. This property is linked to the 
        %   fraction property. Altering one will also alter the other. 
        %   Values must  be finite, non-empty, larger or equal to 0 and 
        %   smaller or equal to 1.        
        %
        %   see also: acoustics.GrainSizeDistribution, fraction
        cumulative_fraction(:,1) double {mustBeFinite, mustBeNonempty, mustBeGreaterThanOrEqual(cumulative_fraction,0) mustBeLessThanOrEqual(cumulative_fraction, 1)}
    end
    methods
        function obj=GrainSizeDistribution(diameter, fraction)
            if nargin>0
                obj.diameter=diameter;
                if nargin>1
                    if numel(diameter)==numel(fraction) % must be cumulative
                        obj.cumulative_fraction=fraction;
                    else
                        obj.fraction=fraction; % must be fraction
                    end
                end
            end
        end
        function set.cumulative_fraction(obj,val)
            obj.fraction=diff(val);
        end
        function val=get.cumulative_fraction(obj)
            val=cumsum([0;obj.fraction]);
        end
        function set.radius(obj,val)
            obj.diameter=val*2;
        end
        function val=get.radius(obj)
            val=obj.diameter.*.5;
        end
        function val=D(obj,N)
            % acoustics.GrainSizeDistribution/D
            % compute percentile of GSD
            %
            %   p=D(obj,N) compute the Nth percentile of the distribution.
            %
            % see also: acoustics.GrainSizeDistribution
            if any(N(:)>1)
                N=N/100;
            end
            obj.check_diameter_fraction();
            [fract,u_id]=unique(obj.cumulative_fraction);
            val=interp1(fract,obj.diameter(u_id),N);
        end
        function val=gsd_int(obj,func,is_radius)
            % acoustics.GrainSizeDistribution/gsd_int
            %   integrates a function of grainsize over the GSD
            %
            %   val=psd_int(obj,func) integrates the function func over the
            %   grain size distribution. If func is not specified, the
            %   function is the grain size, i.e. the mean will be computed
            %
            %   val=psd_int(obj,func,is_radius) is_radius specifies to
            %   integrate a function of the radius instead of the diameter.
            %
            % see also: acoustics.GrainSizeDistribution            
            if nargin < 2 || isempty(func)
                func=@(x) x;
            end
            if nargin < 3 || isempty(is_radius)
                is_radius=false;
            end
            if is_radius
                mult=.5;
            else
                mult=1;
            end
            obj.check_diameter_fraction();
            gs=[obj.diameter];
            val=nansum((func(gs(1:end-1,:)*mult)+func(gs(2:end,:)*mult))*.5.*[obj.fraction],1);
        end
        function val=moment(obj,n,is_radius)
            % acoustics.GrainSizeDistribution/moment
            %   Compute n-th statistical moment of GSD 
            %
            %   val=moment_diam(obj,n) computes the n-th statistical moment
            %   of the GSD. 
            %
            %   val=moment_diam(obj,n,is_radius) computes the n-th
            %   statistical moment of the radius PDF.
            %
            % see also: acoustics.GrainSizeDistribution
            if nargin < 2 || isempty(n)
                n=1;
            end
            if nargin < 3 || isempty(is_radius)
                is_radius=false;
            end
            func=@(x) x.^n;
            val=obj.gsd_int(func,is_radius);
        end
        function val=mean_diam(obj)
            % acoustics.GrainSizeDistribution/mean_diam
            %   Compute the mean diameter
            % 
            %   val=mean_diam(obj) computes the mean diameter of the GSD
            %
            % see also: acoustics.GrainSizeDistribution
            val=obj.moment(1);
        end
        function val=mean_radius(obj)
            % acoustics.GrainSizeDistribution/meam_radius
            %   Compute the mean radius
            %
            %   val=mean_radius(obj) compute the mean radius from the GSD
            %
            % see also: acoustics.GrainSizeDistribution
            val=obj.moment(1,true);
        end
        function ff=form_function(obj, wave_number) 
            % acoustics.GrainSizeDistribution/form_function
            %   Compute the form function of the sediment
            %
            %   ff=form_function(obj, wave_number) computes the form
            %   function or relative backscattering length (Medwin and
            %   Clay) given the GSD and the wave_number. wave_number can be
            %   either the wave number of the acostic signal (1/m) or an
            %   acoustics.PistonTransducer object. This function implements
            %   equation 7 from Thorne and Meral (2012). 
            %
            %   see also: acoustics.GrainSizeDistribution, form_function_static
            narginchk(2,3);
            if isa(wave_number,'acoustics.PistonTransducer')
                wave_number=wave_number.wavenumber;
            end
            ff=acoustics.GrainSizeDistribution.form_function_static(wave_number.*[obj.radius]);
        end
        function chi=norm_scatter_xsection(obj,wave_number)
            % acoustics.GrainSizeDistribution/norm_scatter_xsection
            %   Compute the normalized total scattering cross-section
            %
            %   chi=norm_scatter_xsection(obj, wave_number) computes the 
            %   normalized total scattering cross-section given the GSD and
            %   the wave_number. wave_number can be either the wave number 
            %   of the acostic signal (1/m) or an 
            %   acooustics.PistonTransducer object. This function 
            %   implements equation 9 from Thorne and Meral (2012).
            %
            % see also: acoustics.GrainSizeDistribution
            narginchk(2,3);
            if isa(wave_number,'acoustics.PistonTransducer')
                wave_number=wave_number.wavenumber;
            end
            chi=acoustics.GrainSizeDistribution.norm_scatter_xsection(wave_number.*[obj.radius]);
        end
        function sigbs=bs_xsection(obj,wave_number)
            % acoustics.GrainSizeDistribution/bs_xsection
            %   Compute the mean backscattering cross section
            %
            %   sigbs=bs_xsection(obj,wave_number) computes the mean
            %   backscattering cross-section by integrating the
            %   form-function over the GSD according to equation 4 in Sassi
            %   et al.(2012). wave_number can be either the wave number 
            %   of the acostic signal (1/m) or an 
            %   acooustics.PistonTransducer object.
            %
            % see also: acoustics.GrainSizeDistribution
            if isa(wave_number,'acoustics.PistonTransducer')
                wave_number=wave_number.wavenumber;
            end
            func=@(x) x.^2.*acoustics.GrainSizeDistribution.form_function_static(wave_number*x).^2;
            sigbs=.25*obj.gsd_int(func,true);
        end
        function ks2=ks_squared(obj,wave_number,dens) 
            % acoustics.GrainSizeDistribution/ks_squared
            %  Compute the ks_squared for volume backscattering computation
            %
            %   ks2=ks_squared(obj,wave_number,dens) computes the
            %   ks_squared which can be used to estimate the volume
            %   backscatter strength (Sv=10log10(ks2*Ms)). dens is the
            %   density of the backscattered particles and wave_number can 
            %   be either the wave number of the acostic signal (1/m) or an 
            %   acooustics.PistonTransducer object. Computations are done 
            %   accordint to Sassi et al (2012), eq. 6.
            %
            % see also: acoustics.GrainSizeDistribution
            if isa(wave_number,'acoustics.PistonTransducer')
                wave_number=wave_number.wavenumber;
            end
            dens=reshape(dens,1,[]);
            ks2=obj.bs_xsection(wave_number)*3/4/pi./dens./obj.moment(3,true);
        end
        function zet=specific_atten(obj,wave_number,dens)
            % acoustics.GrainSizeDistribution/specific_atten
            %   Compute the specific attenuation
            %
            %   zet=specific_atten(obj,wave_number,dens) computes the
            %   specific attenuation by integrating the total scattering 
            %   cross-section follwing Sassi et al. 2012, eq 9. wave_number 
            %   can be either the wave number of the acostic signal (1/m) 
            %   or an acooustics.PistonTransducer object. dens is the
            %   density of the scattered particles.
            %
            % see also: acoustics.GrainSizeDistribution
            if isa(wave_number,'acoustics.PistonTransducer')
                wave_number=wave_number.wavenumber;
            end
            dens=reshape(dens,1,[]);
            func=@(x) x.^2.*acoustics.GrainSizeDistribution.norm_scatter_xsection_static(wave_number*x);
            zet=3/4./dens.obj.gsd_int(func,true)./obj.moment(3,true);
        end
        function check_diameter_fraction(obj)
            % acoustics.GrainSizeDistribution/check_diameter_fraction 
            %   Check the diameter and fraction properties
            %
            %   check_diameter_fraction(obj) generates an error if
            %   properties are not ok
            %
            % see also: acoustics.GrainSizeDistribution
            diams=[obj.diameter];
            fract=[obj.fraction];
            assert(size(diams,1)==size(fract,1)+1,'grainsize should have one more element than fraction')           
            assert(all((sum(fract,1)-1).^2 < eps),'sum of fractions in distributions should equal 1')
            assert(all(issorted(diams)),'diameters must be in ascending order')
        end
        function hp=plot_cumulative(obj,varargin)
            % acoustics.GrainSizeDistribution/plot_cumulative
            %   Plots the cumulative GSD
            %
            %   plot_cumulative(obj) creates a plot with the cumulative 
            %   distribution function
            %
            %   plot_cumulative(obj,...) pass options to the plot function
            %
            %   h=plot(...) get handle of plot
            %
            % see also: acoustics.GrainSizeDistribution
            h=plot([obj.diameter],[obj.cumulative_fraction]*100,varargin{:});
            set(gca,'xscale','log')
            xlabel('Diameter (m)')
            ylabel('Cumulative fraction (%)')
            if nargout>0
                hp=h;
            end
        end
        function hp=plot(obj,varargin)
            % acoustics.GrainSizeDistribution/plot
            %   Plots the GSD   
            %
            %   plot_psd(obj) plots the GSD. If obj is not
            %   scalar it will plot the GSDs as lines.
            %
            %   plot_psd(obj,...) pass options for the builtin plot func.
            %
            %   h=plot_psd(obj) return a handle to the plot
            %
            % see also: acoustics.GrainSizeDistribution
            gs=[obj.diameter];
            h=plot(.5*(gs(1:end-1,:)+gs(2:end,:)),[obj.fraction]*100,varargin{:});
            set(gca,'xscale','log')
            xlabel('Diameter (m)')
            ylabel('Fraction (%)')
            if nargout>0
                hp=h;
            end
        end
        function hp=plot_hist(obj,varargin)
            % acoustics.GrainSizeDistribution/plot_hist
            %   Plot a histogram of the GSD
            %   
            %   plot_hist(obj) plot a histogram of the GSD function. 
            %
            %   plot_hist(obj,...) pass options to the builtin patch func.
            %
            %   h=plot_hist(obj) return a handle to the patches
            %
            % see also: acoustics.GrainSizeDistribution   
            if isscalar(obj)
                h=patch(obj.diameter((0:(numel(obj.diameter)-2))+[1;2;2;1;1]), obj.fraction'.*([0;0;1;1;0])*100,zeros(5,numel(obj.fraction)),varargin{:});
            else 
                error('Can only make a histogram of a scalar object')
            end
            set(gca,'xscale','log')
            xlabel('Diameter (m)')
            ylabel('Fraction (%)')
            if nargout>0
                hp=h;
            end
        end
        function plot_form_function(obj,wave_number,varargin)
            % acoustics.GrainSizeDistribution/plot_form_function
            %   plot the form function
            %
            %   plot_form_function(obj,wave_number) plot the form function given the
            %   wave_number.
            %
            %   plot_form_function(obj,wave_number,...) pass options to the builtin plot
            %   function
            %
            %   h=plot_form_function(...) return a handle to the plot
            %
            %   see also: acoustics.GrainSizeDistribution, form_function
            ff=obj.form_function(wave_number,varargin{:});
            gs=[obj.diameter];
            plot(gs,ff)
            set(gca,'xscale','log','yscale','log')
            xlabel('Diameter (m)')
            ylabel('f (-)')
        end
        function plot_all(obj,wave_number,varargin)
            % acoustics.GrainSizeDistribution/plot_all
            %   plot a GSD histogram and the form function
            %
            %   plot_all(obj,wave_number) plot the form function given the
            %   wave_number and a histogram of the GSD.
            %
            %   plot_all(obj,wave_number,...) pass options to the builtin plot
            %   function
            %
            %   h=plot_all(...) return a handle to the plot
            %
            %   see also:  acoustics.GrainSizeDistribution, form_function
            subplot(2,1,1)
            obj.plot_hist;
            subplot(2,1,2)
            obj.plot_form_function(wave_number,varargin{:});
        end
    end
    methods (Static)
        function ff=form_function_static(x) 
            % acoustics.GrainSizeDistribution/form_function_static
            %   computes the form function
            %
            %   acoustics.GrainSizeDistribution.form_function_static(x)
            %   computes the form function given x which is k*a, with k the
            %   acoustic wavenumber and a the particle radius.
            %
            %   see also: acoustics.GrainSizeDistribution, form_function
            ff=x.^2.*(1-0.35*exp(-((x-1.5)/0.7).^2)).*(1+0.5*exp(-((x-1.8)/2.2).^2))./(1+0.9*x.^2);
        end
        function chi=norm_scatter_xsection_static(x)
            % acoustics.GrainSizeDistribution/norm_scatter_xsection_static
            %   computes the normalized scattering cross-section
            %
            %   acoustics.GrainSizeDistribution.form_function_static(x)
            %   computes the form function given x which is k*a, with k the
            %   acoustic wavenumber and a the particle radius.
            %
            %   see also: acoustics.GrainSizeDistribution, norm_scat_xsection
            chi = 0.29*x.^4./(0.95+1.28*x.^2+0.25*x.^4);        
        end
    end
end