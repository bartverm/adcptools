classdef GrainSizeDistribution < handle
    properties
        diameter(:,1) double {mustBeNonnegative, mustBeNonNan, mustBeNonempty} = [0; inf]
        fraction(:,1) double {mustBeFinite, mustBeNonempty, mustBeGreaterThanOrEqual(fraction,0), mustBeLessThanOrEqual(fraction, 1)}=1
    end
    properties(Dependent)
        radius(:,1) double {mustBeNonnegative, mustBeNonNan, mustBeNonempty}
        cumulative_fraction(:,1) double {mustBeFinite, mustBeNonempty, mustBeGreaterThanOrEqual(cumulative_fraction,0) mustBeLessThanOrEqual(cumulative_fraction, 1)}
    end
    methods
        function obj=GrainSizeDistribution(grain_size, fraction)
            if nargin>0
                obj.diameter=grain_size;
                if nargin>1
                    if numel(grain_size)==numel(fraction) % must be cumulative
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
            if any(N(:)>1)
                N=N/100;
            end
            obj.check_grainsize_fraction();
            [fract,u_id]=unique(obj.cumulative_fraction);
            val=interp1(fract,obj.diameter(u_id),N);
        end
        function val=gsd_int(obj,func,is_radius)
            % integrates a function of grainsize over the GSD
            %
            %   val=psd_int(obj,func) integrates the function func over the
            %   grain size distribution. If func is not specified, the
            %   function is the grain size, i.e. the mean will be computed
            %
            %   val=psd_int(obj,func,is_radius) is_radius specifies to
            %   integrate a function of the radius instead of the diameter.
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
            obj.check_grainsize_fraction();
            gs=[obj.diameter];
            val=nansum((func(gs(1:end-1,:)*mult)+func(gs(2:end,:)*mult))*.5.*[obj.fraction],1);
        end
        function val=moment(obj,n,is_radius)
            % Compute n-th statistical moment of GSD 
            %
            %   val=moment_diam(obj,n) computes the n-th statistical moment
            %   of the GSD. 
            %
            %   val=moment_diam(obj,n,is_radius) computes the n-th
            %   statistical moment of the radius PDF.
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
            % Compute the mean diameter
            % 
            %   val=mean_diam(obj) computes the mean diameter of the GSD
            val=obj.moment(1);
        end
        function val=mean_radius(obj)
            % Compute the mean radius
            %
            %   val=mean_radius(obj) compute the mean radius from the GSD
            val=obj.moment(1,true);
        end
        function ff=form_function(obj, wave_number) 
            % Compute the form function of the sediment
            %
            %   ff=form_function(obj, wave_number) computes the form
            %   function or relative backscattering length (Medwin and
            %   Clay) given the GSD and the wave_number. wave_number can be
            %   either the wave number of the acostic signal (1/m) or an
            %   acoustics.PistonTransducer object. This function implements
            %   equation 7 from Thorne and Meral (2012). 
            narginchk(2,3);
            if isa(wave_number,'acoustics.PistonTransducer')
                wave_number=wave_number.wavenumber;
            end
            ff=acoustics.GrainSizeDistribution.form_function_static(wave_number.*[obj.radius]);
        end
        function chi=norm_scatter_xsection(obj,wave_number)
            % Compute the normalized total scattering cross-section
            %
            %   chi=norm_scatter_xsection(obj, wave_number) computes the 
            %   normalized total scattering cross-section given the GSD and
            %   the wave_number. wave_number can be either the wave number 
            %   of the acostic signal (1/m) or an 
            %   acooustics.PistonTransducer object. This function 
            %   implements equation 9 from Thorne and Meral (2012). 
            narginchk(2,3);
            if isa(wave_number,'acoustics.PistonTransducer')
                wave_number=wave_number.wavenumber;
            end
            chi=acoustics.GrainSizeDistribution.norm_scatter_xsection(wave_number.*[obj.radius]);
        end
        function sigbs=bs_xsection(obj,wave_number)
            % Compute the mean backscattering cross section
            %
            %   sigbs=bs_xsection(obj,wave_number) computes the mean
            %   backscattering cross-section by integrating the
            %   form-function over the GSD according to equation 4 in Sassi
            %   et al.(2012). wave_number can be either the wave number 
            %   of the acostic signal (1/m) or an 
            %   acooustics.PistonTransducer object.
            if isa(wave_number,'acoustics.PistonTransducer')
                wave_number=wave_number.wavenumber;
            end
            func=@(x) x.^2.*acoustics.GrainSizeDistribution.form_function_static(wave_number*x).^2;
            sigbs=.25*obj.gsd_int(func,true);
        end
        function ks2=ks_squared(obj,wave_number,dens) 
            %  Compute the ks_squared for volume backscattering computation
            %
            %   ks2=ks_squared(obj,wave_number,dens) computes the
            %   ks_squared which can be used to estimate the volume
            %   backscatter strength (Sv=10log10(ks2*Ms)). dens is the
            %   density of the backscattered particles and wave_number can 
            %   be either the wave number of the acostic signal (1/m) or an 
            %   acooustics.PistonTransducer object. Computations are done 
            %   accordint to Sassi et al (2012), eq. 6.
            if isa(wave_number,'acoustics.PistonTransducer')
                wave_number=wave_number.wavenumber;
            end
            dens=reshape(dens,1,[]);
            ks2=obj.bs_xsection(wave_number)*3/4/pi./dens./obj.moment(3,true);
        end
        function zet=specific_atten(obj,wave_number,dens)
            % Compute the specific attenuation
            %
            %   zet=specific_atten(obj,wave_number,dens) computes the
            %   specific attenuation by integrating the total scattering 
            %   cross-section follwing Sassi et al. 2012, eq 9. wave_number 
            %   can be either the wave number of the acostic signal (1/m) 
            %   or an acooustics.PistonTransducer object. dens is the
            %   density of the scattered particles.
            if isa(wave_number,'acoustics.PistonTransducer')
                wave_number=wave_number.wavenumber;
            end
            dens=reshape(dens,1,[]);
            func=@(x) x.^2.*acoustics.GrainSizeDistribution.norm_scatter_xsection_static(wave_number*x);
            zet=3/4./dens.obj.gsd_int(func,true)./obj.moment(3,true);
        end
        function n=n_particles(obj,vol_concentration)
            n=vol_concentration*3/4/pi/obj.moment(3,true);
        end
        function check_grainsize_fraction(obj)
            diams=[obj.diameter];
            fract=[obj.fraction];
            assert(size(diams,1)==size(fract,1)+1,'grainsize should have one more element than fraction')           
            assert(all((sum(fract,1)-1).^2 < eps),'sum of fractions in distributions should equal 1')
            assert(all(issorted(diams)),'diameters must be in ascending order')
        end
        function plot(obj)
            plot([obj.diameter],[obj.cumulative_fraction]*100,'.-')
            set(gca,'xscale','log')
            xlabel('Diameter (m)')
            ylabel('Cumulative fraction (%)')
        end
        function plot_psd(obj)
            if isscalar(obj)
                patch(obj.diameter((0:(numel(obj.diameter)-2))+[1;2;2;1;1]), obj.fraction'.*([0;0;1;1;0])*100,zeros(5,numel(obj.fraction)));
            else
                gs=[obj.diameter];
                plot(.5*(gs(1:end-1,:)+gs(2:end,:)),[obj.fraction]*100);
            end
            set(gca,'xscale','log')
            xlabel('Diameter (m)')
            ylabel('Fraction (%)')
        end
        function plot_form_function(obj,varargin)
            ff=obj.form_function(varargin{:});
            gs=[obj.diameter];
            plot(gs,ff)
            set(gca,'xscale','log','yscale','log')
            xlabel('Diameter (m)')
            ylabel('f (-)')
        end
        function plot_all(obj,varargin)
            subplot(2,1,1)
            obj.plot_psd;
            subplot(2,1,2)
            obj.plot_form_function(varargin{:});
        end
    end
    methods (Static)
        function ff=form_function_static(x) % Thorne and Meral, 2012, eq.7 
            ff=x.^2.*(1-0.35*exp(-((x-1.5)/0.7).^2)).*(1+0.5*exp(-((x-1.8)/2.2).^2))./(1+0.9*x.^2);
        end
        function chi=norm_scat_xsection_static(x) % Thorne and Meral, 2012, eq. 9
            chi = 0.29*x.^4./(0.95+1.28*x.^2+0.25*x.^4);        
        end
    end
end