classdef Sv2SSC_Power< handle
% Compute SSC from ADCP  backscatter using the power-law method
%
%   obj=acoustics.Sv2SSC_Power();
%
%   acoustics.Sv2SSC_Power properties:
%   samples - water samples for calibration
%   sv - backscatter values corresponding to samples for calibration
%
%   acoustics.Sv2SSC_Power methods:
%   calibrate_a_b - calibrate the power law constants
%   mass_concentration - compute mass concentration from backscatter
%
%   see also: acoustics
    properties
        samples(:,1) acoustics.WaterSample
        sv(:,1) double
    end
    methods
        function [a,b]=calibrate_a_b(obj)
            assert(numel(obj.samples)==numel(obj.sv))
            ssc=log10([obj.samples.mass_concentration]/1000); % mass concentration in mg/L
            pars=robustfit(obj.sv,ssc);
            b=pars(1);
            a=pars(2);
        end
        function ssc=mass_concentration(obj,sv)
            [a,b]=obj.calibrate_a_b();
            ssc=1e3*10^(a*sv+b); % output kg/m3
        end
    end
end