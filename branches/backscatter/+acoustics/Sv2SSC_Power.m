classdef Sv2SSC_Power < acoustics.Sv2SSC
% Compute SSC from ADCP  backscatter using the power-law method
%
%   obj=acoustics.Sv2SSC_Power() Construct default object
%
%   obj=acoustics.Sv2SSC_Power(objects) Pass objects to assign properties of
%   calibration. Objects of type ADCP are assigned to the 'adcp' property,
%   objects of type acoustics.WaterSample are assigned to the property
%   'samples' and duration objects are assigned to the averaging period
%   object.
%
%   Sv2SSC_Power properties:
%   adcp - ADCP objects used for calibration
%   samples - acoustics.WaterSample objects for calibration
%   averaging_period - duration object to indicate averaging period
%
%   Sv2SSC_Power methods:
%   mass_concentration - mass concentration in g/L
%   match_adcp_samples - match samples with corresponding ADCP data
%   sv_for_calibration - averaged backscatter strength for calibration
%   plot - plot the backscatter strength vs mass concentration
%   calibrate_a_b - calibrate the power law constants
%
%   see also: ADCP, acoustics, Sv2SSC_Power, Sv2SSC_Sassi
    methods
        function [a,b,stda, stdb]=calibrate_a_b(obj)
        % Calibrates constants of power law method
        %
        %   [a, b, stda, stdb]=calibrate_a_b(obj) Calibrate the constants
        %   a and b for the power law method. stda and stdb hold the
        %   standard deviation of these constants.
        %
        %   see also: Sv2SSC_Power, ADCP, acoustics
            sv_cal=reshape(obj.sv_for_calibration(),[],1);
            ssc=reshape(log10([obj.samples.mass_concentration]/1000),[],1); % mass concentration in mg/L
            [pars, stdpars]=lscov([sv_cal ones(size(sv_cal))],ssc);
            a=pars(1);
            stda=stdpars(1);
            b=pars(2);
            stdb=stdpars(2);
        end
        function ssc=mass_concentration(obj,sv)
        % Computes the mass concentration given backscatter strength
        %
        %   ssc=mass_concentration(obj) computes the suspended sediment
        %   concentration in g/L using the backscatter strength for the
        %   ADCP objects.
        %
        %   ssc=mass_concentration(obj,sv) specify the backscatter strength
        %   to be used for sediment concentration computation
        %
        %   see also: Sv2SSC_Power, Sv2SSC, acoustics, ADCP
            if nargin < 2
                sv=obj.adcp.backscatter();
            end
            [a,b]=obj.calibrate_a_b();
            ssc=1e3*10.^(a*sv+b); % output kg/m3
        end
        function plot(obj)
        % Plot the calibration
        %
        %   Plots the backscatter strength vs the sediment concentration
        %   and adds a line representing the calibrated parameters.
        %
        %   see also: Sv2SSC_Power, Sv2SSC, acoustics, ADCP
            plot@acoustics.Sv2SSC(obj);
            hold on
            x=get(gca,'xlim');
            y=obj.mass_concentration(x);
            plot(x,y,'r-')
        end
    end
end