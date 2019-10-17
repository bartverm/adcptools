classdef GrainSizeDistribution < handle
    properties
        grainsize(:,1) double {mustBeNonnegative, mustBeNonNan, mustBeNonempty} = [0; inf]
        fraction(:,1) double {mustBeFinite, mustBeNonempty, mustBeGreaterThanOrEqual(fraction,0), mustBeLessThanOrEqual(fraction, 1)}=1
    end
    properties(Dependent)
        cumulative_fraction(:,1) double {mustBeFinite, mustBeNonempty, mustBeGreaterThanOrEqual(cumulative_fraction,0) mustBeLessThanOrEqual(cumulative_fraction, 1)}
    end
    methods
        function obj=GrainSizeDistribution(grain_size, fraction)
            if nargin>0
                obj.grainsize=grain_size;
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
        function val=D(obj,N)
            if any(N(:)>1)
                N=N/100;
            end
            obj.check_grainsize_fraction();
            [fract,u_id]=unique(obj.cumulative_fraction);
            val=interp1(fract,obj.grainsize(u_id),N);
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
            val=nansum((func(obj.grainsize(1:end-1)*mult)+func(obj.grainsize(2:end)*mult))*.5.*obj.fraction);
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
            val=obj.moment(1);
        end
        function val=mean_radius(obj)
            val=obj.moment(1,true);
        end
        function [ff, x]=form_function(obj, wave_number, radius) % Thorne and Meral, 2008
            narginchk(2,3);
            if isa(wave_number,'acoustics.PistonTransducer')
                wave_number=wave_number.wavenumber;
            end
            if nargin < 3
                radius=obj.grainsize/2;
            end
            x=wave_number.*radius;
            ff=x.^2.*(1-0.35*exp(-((x-1.5)/0.7).^2)*1+0.5*exp(-((x-1.8)/2.2).^2))./(1+0.9*x.^2);
        end
        function sigbs=bs_xsection(obj,wave_number) % Medwin and Clay, 1998
            func=@(x) x.^2.*obj.form_function(wave_number,x).^2;
            sigbs=.25*obj.gsd_int(func,true);
        end
        function n=n_particles(obj,vol_concentration)
            n=vol_concentration*3/4/pi/obj.moment(3,true);
        end
        function check_grainsize_fraction(obj)
            assert(numel(obj.grainsize)==numel(obj.fraction)+1,'grainsize should have one more element than fraction')           
            assert(sum(obj.fraction)-1 < eps,'sum of all fractions should equal 1')
            assert(issorted(obj.grainsize),'grainsize must be in ascending order')
        end
        function plot(obj)
            plot(obj.grainsize,obj.cumulative_fraction*100,'.-')
            set(gca,'xscale','log')
            xlabel('Diameter (m)')
            ylabel('Cumulative fraction (%)')
        end
        function plot_psd(obj)
            patch(obj.grainsize((0:(numel(obj.grainsize)-2))+[1;2;2;1;1]), obj.fraction'.*([0;0;1;1;0])*100,zeros(5,numel(obj.fraction)));
            set(gca,'xscale','log')
            xlabel('Diameter (m)')
            ylabel('Fraction (%)')
        end
        function plot_form_function(obj,varargin)
            ff=obj.form_function(varargin{:});
            gs=obj.grainsize();
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
end