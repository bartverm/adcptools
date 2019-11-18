classdef Sv2SSC_ConstantGSD < acoustics.Sv2SSC
% Sv2SSC_ConstantGSD computes mass conc. from backscatter (Sassi et al., 2012)
%
%   obj=Sv2SSC_ConstantGSD() Constructs default object
%
%   Sv2SSC_ConstantGSD properties:
%   adcp - ADCP objects used for calibration
%   samples - acoustics.WaterSample objects for calibration
%   reference_idx - index to select reference samples
%   attenuation_idx - index to select attenuation samples 
%   averaging_period - duration object to indicate averaging period
%
%   Sv2SSC_ConstantGSD methods:
%   mass_concentration - mass concentration in g/L (abstract method)
%   match_adcp_samples - match samples with corresponding ADCP data
%   sv_for_calibration - averaged backscatter strength for calibration
%   plot - plot the backscatter strength vs mass concentration
%
%   Sv2SSC_ConstantGSD methods:
%   calibrate_b_gamma - calibrate the b and gamma calibration constants
%   mass_concentration - estimate the mass concentration profiles
%
%   Sv2SSC_ConstantGSD static methods:
%   sv_to_beta - compute beta from Sv values
%   delta_Sv - find correction for Sv based on new Kc and C
%
%   see also: ADCP, acoustics, Sv2SSC_Power, Sv2SSC

    properties
        % Sv2SSC_ConstantGSD/reference_idx
        %
        % indices pointing to the water samples to be used as reference
        % samples, i.e. collected near the tranducer head
        %
        % see also: Sv2SSC_ConstantGSD, acoustics, ADCP
        reference_idx(:,1) double {mustBeFinite, mustBeInteger, mustBePositive} = double.empty(0,1);
        
        % Sv2SSC_ConstantGSD/attenuation_idx
        %
        % indices pointing to the water samples to be used for attenuation
        % estimates
        %
        % see also: Sv2SSC_ConstantGSD, acoustics, ADCP
        attenuation_idx(:,1) double {mustBeFinite, mustBeInteger, mustBePositive} = double.empty(0,1);      
    end
    methods        
%         %%% Ordinary methods %%%
        function [b, stdb, gamma_e, xi_e]=calibrate_b_gamma(obj)
        % Estimate the calibration parameters b and gamma
        %
        %   [b, stdb, gamma_e, xi_e]=calibrate_b_gamma(obj) calibrates the
        %   parameters b and gamma. For b the reference samples are used
        %   and a single value is returned. stdb holds the standard
        %   deviation of the estimated parameter. gamma_e is based on the
        %   attenuation samples. For each sample a value of gamma_e is
        %   returned. xi_e is computed from gamme_e and can be used to
        %   estimated attenuation coefficients.
        %
        %   see also: Sv2SSC_ConstantGSD, acoustics
        
            sv=obj.sv_for_calibration();
            ssc=[obj.samples.mass_concentration];
            [~, cell_idx, ensemble_idx]=obj.match_adcp_samples();
            srange=obj.adcp.depth_cell_slant_range;
            sv_all=obj.adcp.backscatter;
            beta_all=acoustics.Sv2SSC_ConstantGSD.sv_to_beta(sv_all);
            % calibration of b, for reference backscatter (paragraph 14)
            
            Sv_ref=reshape(sv(obj.reference_idx),[],1);
            beta_ref=acoustics.Sv2SSC_ConstantGSD.sv_to_beta(Sv_ref); % compute beta for ref (eqs. 12, Sassi et al.)
            ms_ref=reshape(ssc(obj.reference_idx),[],1);
            k_ref=beta_ref./ms_ref; % compute k for ref (eqs. 12, Sassi et al.)
            [b,stdb]=lscov(Sv_ref,10*log10(k_ref)); % simple linear fitting forced through the origin
            
            % calibration of specific attenuation (paragraph 15)
            Sv_cal=reshape(sv(obj.attenuation_idx),[],1);
            beta_cal=acoustics.Sv2SSC_ConstantGSD.sv_to_beta(Sv_cal);
            ms_cal=reshape(ssc(obj.attenuation_idx),[],1);
            [gamma_e,k_ref,int_beta]=deal(nan(numel(obj.attenuation_idx),1));
            for cs=1:numel(obj.attenuation_idx)
                ensid=ensemble_idx{obj.attenuation_idx(cs)};
                k_ref(cs)=nanmean(beta_all(1,ensid,:),[2,3]).^b;
                beta_cur=nanmean(beta_all(:,ensid,:),[2,3]);
                srange_cur=nanmean(srange(:,ensid,:),[2,3]);
                cell_idx_cur=mode(cell_idx{obj.attenuation_idx(cs)},[2,3]);
                int_beta(cs)=sum((beta_cur(1:cell_idx_cur-1)+beta_cur(2:cell_idx_cur)).*diff(srange_cur(1:cell_idx_cur))/2); % trapezoidal integration
                gamma_e(cs)=(k_ref(cs)-beta_cal(cs)/ms_cal(cs))./int_beta(cs);
            end
            xi_e=gamma_e*5/log(10);
        end
        function val=mass_concentration(obj,inadcp,b,gamma_e)
        % Estimate ssc based on the calibration
        %
        %   ssc=mass_concentration(obj) Estimates the mass concentration in
        %   g/L (kg/m^3) for the backscatter from obj.adcp. 
        %
        %   ssc=mass_concentration(obj,inadcp) allows to specify the adcp
        %   object for which the mass concentration is to be computed
        %
        %   ssc=mass_concentration(obj,inadcp, b, gamma_e) to additionally
        %   specify values of b and gamma. Default is to use the output
        %   from the calibrate_b_gamma method, where the value for gamma is
        %   taken to be the median of the estimated values.
        %
        %   see also: Sv2SSC_ConstantGSD, acoustics
        
            if nargin < 4
                [b_tmp,gamma_e_tmp]=obj.calibrate_b_gamma();
                gamma_e=nanmedian(gamma_e_tmp);
                if nargin < 3
                    b=b_tmp;
                    if nargin<2
                        inadcp=obj.adcp;
                    end
                end
            end
            Sv=inadcp.backscatter();
            srange=inadcp.depth_cell_slant_range();
            beta=acoustics.Sv2SSC_ConstantGSD.sv_to_beta(Sv);
            K_ref=beta(1,:,:).^b; % Kref from calibration

            int_beta=[zeros(1,size(beta,2), size(beta,3)); cumsum((beta(1:end-1,:,:)+beta(2:end,:,:)).*diff(srange,1,1)/2,1)];
            val=beta./(K_ref-gamma_e.*int_beta);
        end
        function plot(obj)
        % make a plot for the calibration of b
        %
        % see also: Sv2SSC_ConstantGSD, acoustics
            sv=obj.sv_for_calibration();
            ssc=[obj.samples.mass_concentration];
            Sv_ref=reshape(sv(obj.reference_idx),[],1);
            beta_ref=acoustics.Sv2SSC_ConstantGSD.sv_to_beta(Sv_ref); % compute beta for ref (eqs. 12, Sassi et al.)
            ms_ref=reshape(ssc(obj.reference_idx),[],1);
            k_ref=beta_ref./ms_ref;
            b=lscov(Sv_ref,10*log10(k_ref)); % simple linear fitting forced through the origin
            plot(Sv_ref, k_ref,'o');
            set(gca,'yscale','log')
            xlabel('S_{v,\alpha_s=0} (dB)')
            ylabel('K(R_{ref}) (-)')
            hold on
            xl=get(gca,'xlim');
            plot(xl,acoustics.Sv2SSC_ConstantGSD.sv_to_beta(xl).^b,'r')
        end
    end
    methods(Static)
        function beta=sv_to_beta(sv)
        % Estimate beta from backscatter values (equation 12a)
            beta=10.^(sv/10);
        end
    end
    
end


