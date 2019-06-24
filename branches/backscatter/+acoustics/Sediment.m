classdef Sediment < handle
% Sediment class
%   Query objects of this class to obtain percentiles of the grain size 
%   distribution, the parameters of the normal phi (-log2(diameter in mm))
%   distribution and the parameters of the lognormal distribution of the
%   sediment diameters in m.
%
%   OBJ=Sediment Constructs a default object with the properties set to
%           their default values
%   OBJ=Sediment(DISTRIB) Creates an object with the given distribution.
%           For the format of this variable see Properties
%   OBJ=Sediment(DISTRIB,DENSITY) Creates and object with given
%           distribution and density
%
%   Sediment Properties:
%       density     - Density of specific weight of the sediment in kg/m^3
%                     or g/L. Defualt is 2650
%       distrib     - m by 2 matrix defining the known percentiles of the 
%                     sediment distribution. The first column contains the
%                     percantage of sediment with a size smaller than the
%                     size given in the second collumn in m. The matrix
%                     must contain at least two rows. The default is: [50
%                     300e-6; 90 450e-6]. This means: D50 = 300 micron,
%                     D90=450 micron. 
%       normParms   - Paramters (mu and sigma) of the normal distribution 
%                     in phi units (Get only property)
%       lognParms   - Parameters (mu and sigma) of the lognormal 
%                     distribution in m (get only property)
%
%   Sediment Methods:
%       plot        - Plots the fitted sediment distribution and the given
%                     quantiles in DISTRIB
%       D(n)        - Returns the size in m for a certain quantile. n in
%                     percentage
%       lognM(N)    - Returns the Nth moment of the lognormal distribution
%

%   Author: Bart Vermeulen
%   Last edit: 06_2019

    properties
        density=2650 %Sediment density in kg/m^3
    end
    properties(SetObservable, AbortSet)
        distrib=[50 300e-6;...  % Grain size distribution (quantiles in %, second column size in m) 
                 90 450e-6]
    end
    properties(Dependent)
        normParms  % Parameters for the normal phi probability distribution
        E
        VAR
    end
    properties(GetAccess=public, SetAccess=private)
        lognParms   % Parameters for the lognormal grain size distribution in m
    end
    methods
        %% Constructor
        function obj=Sediment(varargin)
            if nargin == 0
                setgrdist(obj);
            end
            addlistener(obj,'distrib','PostSet',@obj.setgrdist);
            if nargin > 0, obj.distrib=varargin{1};end
            if nargin == 2,obj.density=varargin{2};end
            assert(nargin<3,'Too many input arguments');
        end
        %% Propety access methods
        function value=get.normParms(obj)
            value(1)=-((obj.lognParms(1)+3*log(10))/log(2)+1);
            value(2)=obj.lognParms(2)/log(2);
        end
        function value=get.E(obj)
            value=obj.lognM(1);
        end
        function value=get.VAR(obj)
            value=obj.lognM(2);
        end
        function set.density(obj,value)
            if isnumeric(value) && isscalar(value) && value>0
                obj.density=value;
            else
                error('sediment:setDensity:BadInput','Density should be numerical positive scalar')
            end
        end
        function set.distrib(obj,value)
            if isnumeric(value) && size(value,2)==2 && size(value,1)>1 && issorted(value(:,1)) && issorted(value(:,2)) && all(value(:)>0) && all(value(:,1)<100)
                obj.distrib=value;
            else
                error('sediment:setDistrib:BadInput','Wrong format for distrib')
            end
        end
        %% Generic methods
        function value=lognM(obj,n)
            % OBJ.lognM(N) Determines the Nth central moment of the lognormal
            %   distribution of grain size distribution in m
            
            % Last edit: 14-04-2010
            mu=obj.lognParms(1);
            sigma=obj.lognParms(2);
            switch n
                case 1
                    value=exp(mu+sigma^2/2);
                case 2
                    value=(exp(sigma^2)-1)*exp(2*mu+sigma^2);
                case 3
                    value=(exp(sigma^2)+2)*sqrt(exp(sigma^2)-1);
                case 4
                    value=exp(4*sigma^2)+2*exp(3*sigma^2)+3*exp(2*sigma^2)-3;
                otherwise
                    value=nan;
            end
        end
        function value=lognrawM(obj,n)
            % OBJ.lognrawM(N) Determines the Nth raw moment of the lognormal
            %   distribution of grain size distribution in m
            
            % Last edit: 14-04-2010
            mu=obj.lognParms(1);
            sigma=obj.lognParms(2);
            value=exp(n*mu+n^2*sigma^2/2);
        end

        function value=D(obj,n) 
            % OBJ.D determines quantiles of the sediment grain size
            % distribution
            %
            % SIZE=OBJ.D(N) Returns the sediment size such that N % of the 
            %       sediment is smaller. E.g. OBJ.D(50) returns the median
            %       grain size for the estimated distribution
            %
          
            % Last edit 13-04-2010
            value=logninv(n/100,obj.lognParms(1),obj.lognParms(2));
        end
        function plot(obj)
           % OBJ.PLOT plots the estimated cdf and the given distribution
           %
           % OBJ.plot Creates the plot
           
           % Last edit: 13-04-2010
           x=linspace(obj.D(0.01),obj.D(99.9),1000);
           y=logncdf(x,obj.lognParms(1),obj.lognParms(2));
           plot(-log2(x*1000),y)
           hold on
           plot(-log2(obj.distrib(:,2)*1000), obj.distrib(:,1)/100,'rx');
           set(gca,'xdir','reverse')
           xlabel('Grain size in phi units \leftarrow')
           ylabel('Cumulative probability \rightarrow')
           hold off
        end
    end
    %% Private methods
    methods(Access=private)
        function perf=lgnfit(obj,x)
            % Determines sum of squared differences between known quantiles
            % and the fit
            if x(2)<0 
                perf=Inf;
            else
                perf=nansum((obj.distrib(:,1)/100-logncdf(obj.distrib(:,2),x(1),x(2))).^2);
            end
        end
        function setgrdist(obj,~,~)
            % Estimates the parameters of the lognormal grain size
            % distribution
            D84=interp1(obj.distrib(:,1),obj.distrib(:,2),84.1,'pchip','extrap')*1000;
            D16=interp1(obj.distrib(:,1),obj.distrib(:,2),15.9,'pchip','extrap')*1000;
            muphi=-log2(sqrt(D84*D16));
            sigphi=log2(sqrt(D84/D16));
            x0(1)=-muphi*log(2)-3*log(10);
            x0(2)=sigphi*log(2);
            [obj.lognParms,~,xflag,outpt]=fminsearch(@obj.lgnfit,x0,optimset('Display','off','TolX',1e-6));
            if xflag ~=1
                warning('sediment:setgrdist:BadEstimate',['Could not find a good distribution fit: ',outpt.message])
            end
        end
    end
    
end