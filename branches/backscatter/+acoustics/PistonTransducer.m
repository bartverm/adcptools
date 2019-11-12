classdef PistonTransducer < handle
    % acoustics.PistonTransducer to get properties of piston transducer.
    %   This class aids in computing characteristics of piston transducers.
    %   Most computations are based on the book by Medwin and Clay.
    % 
    % acoustics.PistonTransducer properties:
    %   radius - radius of transducer in m (default is 0.0505)
    %   frequency - frequency of the transducer in Hz (default is 614.4e3)
    %   depth - depth of transducer in m (default is 1)
    %   water - acoustics.Water object describing properties of water
    %   beam_width - 3dB beam width in radians 
    %
    % acoustics.PistonTransducer dependent properties (read only):
    %   wavelength - wavelength in m
    %   wavenumber - angular wave number of the sound waves in 1/m
    %   nearfield - range of near field in m
    %   angularfreq - angular frequency of the sound waves in rad/s
    %   attenuation - sound attenuation caused by water in dB/m
    %   attenuation_e - sound attenuation caused by water in 1/m (neper)
    %   speedsound - speed ouf sound in water in m/s
    %
    % acoustics.PistonTransducer methods:
    %   directivity - Directivity of the emitted sound (neper)
    %   dir_response - Directional response (dB) of the transducer 
    %   near_field_correction - Compute the near field correction factor
    %   side_lobes - Compute angle of transmitted side lobes
    %   zeros - Compute angle (radians) where no sound it transmitted
    %   plot - Plot the directional response of the transducer
    %   plot_3d - Make a 3D plot of the directional response
    
    properties
        % radius - Radius of the transducer in m
        %   This should be a positive, double, Default: 0.0505
        %
        radius double = 0.0505 
        
        % frequency - Frequency of the emitted sound in Hz
        %   This should be a positive, double, scalar value. Default:
        %   614.4e3
        frequency double = 614.4e3 
        
        % depth - Depth of the transducer in m
        %   This should be a positive, double, scalar value. Default: 0
        depth double = 0
        
        % water - Properties of water
        %   This is a scalar object of class acoustics.Water. Holds common
        %   properties of water
        %
        %   see also: acoustics.Water
        water(1,1) acoustics.Water
    end
    properties(Dependent, GetAccess=public, SetAccess=private)
        wavelength %The wave length of the sound produced by the transducer in m.
        wavenumber %The (angular) wave number of the sound produced by the transducer in 1/m.
        nearfield %The near field range in m. After this range the planar wave assumption holds.
        angularfreq % The angular frequency of the sound produced by the transducer in rad/s.
        attenuation % Sound attenuation in water in dB/m
        attenuation_e % Sound attenuation in 1/m (neper)
        speedsound % Speed of sound in the given water m/s
        beam_width % -3dB beam width
    end
    properties(Dependent, Access=private)
       f1,f2,A1,A2,A3,P2,P3
    end
    properties(Constant, Access=private)
        P1=1;                                                               % -, Medwin and clay, pg 110, eq 3.4.30
        besjzeros=[3.83171 7.01559 10.17347 13.32369 16.47063 19.61586...
            22.76008 25.90367 29.04683 32.18968 35.33231 38.47477...
            41.61709 44.75932 47.90146 51.04354 54.18555 57.32753...
            60.46946 63.61136];                                             % Zeros of bessel function j1 (from Abramowitz and Stegun table 9.5)
    end
    methods(Access=public)
        function Di=directivity(obj,phi)
            % Directivity of the transducer 
            % 
            % DI=obj.directivity(PHI) computes the directivity (-) as a function 
            %   of the angle phi  (radians) from the transducer axis. Output
            %   DI will have same size as input PHI.
            assert(isnumeric(phi),'Phi should be numerical')
            Di=2*besselj(1,obj.wavenumber*obj.radius.*sin(phi))./(obj.wavenumber*obj.radius.*sin(phi));
        end
        function Dr=dir_response(obj,phi)
            % Directional response of the transducer
            %
            % DR=obj.dir_response(PHI) returns the directional responce in 
            %   dB of the transducer as a function of the angle PHI (rad) 
            %   from the transducer axis. Output DR has the same size as
            %   input PHI.
            assert(isnumeric(phi),'Phi should be numerical')
            Dr=10*log10(obj.directivity(phi).^2);
            Dr(phi==0)=0;
        end       
        function psi=near_field_correction(obj,range)
            % Returns the near field correction factor for the transducer
            %
            % PSI=obj.near_field_correction(RANGE) returns the near field
            %  correction factor PSI (-) (Downing et al. 1995), to correct for 
            %  interference in the nearfield of the transducer as a function of
            %  the RANGE from the transducer. Output PSI will have the same
            %  size as input RANGE
            assert(isnumeric(range),'Range should be numerical')
            z=range./obj.nearfield;
            psi=(1+1.32*z+(2.5*z).^3.2)./(1.32*z+(2.5*z).^3.2);
        end
        function ang=zeros(obj,n)
            % Computes the angles at which no sound is transmitted
            %
            % ang=obj.zeros(n) with n being the n-th zero, returns the
            % angle of the n-th zero from the axis of the acosutic beam
            % in radians
            assert(isnumeric(n) && all(uint32(n(:))==n(:)) && all(n(:)>0) && all(n(:)<21),'n must be an integer between 1 and 20')
            ang=asin(obj.besjzeros(n)./obj.wavenumber./obj.radius);
        end
        function ang=side_lobes(obj,n)
            % Computes the angles at which side lobes are transmitted
            %
            % ang=obj.side_lobes(n) with n being the n-th size lobe, returns the
            % angle of the n-th zero from the axis of the acosutic beam
            % in radians. For n=0 return the angle of the main lobe, i.e. 0
            assert(isnumeric(n) && all(uint32(n(:))==n(:)) && all(n(:)>=0) && all(n(:)<20),'n must be an integer between 0 and 19')
            ang=nan(size(n));
            ang(n==0)=0;
            for ci=unique(n(n>0))
                ang(n==ci)=fminbnd(@(x) -(obj.directivity(x).^2),obj.zeros(ci),obj.zeros(ci+1),optimset('TolX',1e-6));
            end
        end
        function plot(obj)
            x=linspace(-pi/2,pi/2,2000);
            hold_stat=get(gca,'NextPlot');
            plot(x/pi*180, obj.dir_response(x))
            hold on
            plot(obj.beam_width*[-1 1]*.5/pi*180, [-3 -3],'r')
            text(obj.beam_width/pi*180,-3,['Beam width: ', num2str(obj.beam_width/pi*180),'^\circ'])
            set(gca,'NextPlot',hold_stat,'Ylim', [-70 0])
            xlabel('Angle (^\circ)')
            ylabel('Directional response (dB)')
        end
        function plot_3d(obj)
            nphi=1000;
            ntht=200;
            db_cutoff=70;
            db_step=10;
            n_lobes=10;
            n_cols=floor(db_cutoff/db_step);
            z=obj.zeros(1:n_lobes);
            phi=linspace(-z(end), z(end),nphi);
            phi=sort([phi, -z(1:end-1), z(1:end-1)]);
            tht=linspace(0,2*pi,ntht);
            [PHI,THETA]=meshgrid(phi,tht);
            DIR=max(obj.dir_response(PHI)+db_cutoff, 0);
            [X,Y,Z]=sph2cart(THETA,pi/2-PHI,DIR);
            surf(X,Y,Z,DIR-db_cutoff)
            shading interp
            lightangle(0,45)
            lighting gouraud
            axis off
            colormap(parula(n_cols))
            set(gca,'clim',[-db_step*n_cols 0],'projection','perspective')
            hc=colorbar('SouthOutside');
            xlabel(hc,'Directional Response (dB)')
            set(hc,'xtick',-db_step*n_cols:db_step:0)
       end
    end
    methods        
        %% Property get methods for public, dependent properties
        function value=get.speedsound(obj)
            value=obj.water.speedsound(obj.depth);
        end
        function value=get.wavelength(obj)
            value=obj.speedsound/obj.frequency;
        end
        function value=get.wavenumber(obj)
            value=2*pi*obj.frequency./obj.speedsound;
        end
        function value=get.nearfield(obj)
            value=pi*obj.radius^2./obj.wavelength;
        end
        function value=get.angularfreq(obj)
            value=2*pi*obj.frequency;
        end
        function value=get.attenuation(obj)
           freq=obj.frequency/1000;
           value=(obj.A1.*obj.P1.*obj.f1.*freq.^2./(freq.^2+obj.f1.^2)+...
                  obj.A2.*obj.P2.*obj.f2.*freq.^2./(freq.^2+obj.f2.^2)+...
                  obj.A3.*obj.P3.*freq.^2)/1000;                     % dB/m, Medwin and clay, pg 109, eq 3.4.29
        end
        function value=get.attenuation_e(obj)
            value=obj.attenuation/8.68;                                    % Medwin and clay, pg104, eq 3.4.8b
        end
        function val=get.beam_width(obj)
            val=2*abs(fzero(@(x) obj.dir_response(x)+3,0));
        end
        %% Property get methods for protected or private properties
        % Boric acid components in Sea Water
        function value=get.A1(obj)
            value=8.68./obj.speedsound.*10.^(0.78*obj.water.pH-5);               % dB/km/kHz, Medwin and clay, pg 110, eq 3.4.30
        end
        function value=get.f1(obj)
            S=obj.water.salinity/1000;
            value=2.8*(S/35).^.5.*10.^(4-1245./(273+...
                obj.water.temperature));                                         % kHz, Medwin and clay, pg 110, eq 3.4.30
        end
        % Magnesium Sulfate components in Sea Water
        function value=get.A2(obj)
            S=obj.water.salinity/1000;
            value=21.44*S./obj.speedsound.*(1+...
                0.025*obj.water.temperature);                                    % dB/km/kHz, Medwin and clay, pg 110, eq 3.4.31
        end
        function value=get.P2(obj)
            value=1-1.37e-4.*obj.depth+6.2e-9.*obj.depth.^2;               % - , Medwin and clay, pg 110, eq 3.4.31
        end   
        function value=get.f2(obj)
            S=obj.water.salinity/1000;
            value=8.17*10.^(8-1990./(273+obj.water.temperature))./(1+...
                0.0018*(S-35));                                 % kHz., Medwin and clay, pg 110, eq 3.4.31
        end
        %Pure Water constants
        function value=get.A3(obj)
            value=obj.speedsound*nan;
            fcold=obj.water.temperature<=20;
            value(fcold)=4.937e-4-2.59e-5*obj.water.temperature(fcold)+...
                9.11e-7*obj.water.temperature(fcold).^2-...
                1.5e-8*obj.water.temperature(fcold).^3;                          % dB/km/kHz, Medwin and clay, pg 110, eq 3.4.32
            fcold=~fcold;
            value(fcold)=3.964e-4-1.146e-5*obj.water.temperature(fcold)+...
                1.45e-7*obj.water.temperature(fcold).^2-...
                6.5e-10*obj.water.temperature(fcold).^3;                         % dB/km/kHz, Medwin and clay, pg 110, eq 3.4.33
        end
        function value=get.P3(obj)
            value=1-3.83e-5*obj.depth+4.9e-10*obj.depth.^2;                % -, Medwin and clay, pg 110, eq 3.4.34
        end
    end
end