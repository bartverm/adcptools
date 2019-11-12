classdef Sv2SSC_ConstantGSD < Sv2SSC
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
%   see also: ADCP, acoustics, Sv2SSC_Power, Sv2SSC_ConstantGSD

% TODO: work with ADCP objects, make Kc_C a dependent property similar to
% sediment_density, testing!
    methods        
%         %%% Ordinary methods %%%
%         function [Kc_C_out, Sv_new]=Kc_C(obj,Sv_cal,Echo_cal,distance_cal,Sv, Echo)
%             if obj.calibrate_Kc_C % if we want the calibration
% 
%                 % check number of samples dimensions match
%                 samps=obj.reference_samples;
%                 nsamps=numel(samps);
%                 assert(size(Sv_cal,2)==nsamps, ['Number of reference samples (',num2str(nsamps),') should match with number of columns (',num2str(size(Sv_cal,2)),') for Sv_cal input'])
%                 assert(size(Echo_cal,2)==nsamps, ['Number of reference samples (',num2str(nsamps),') should match with number of columns (',num2str(size(Echo_cal,2)),') for Echo_cal input'])
%                 assert(size(distance_cal,2)==nsamps, ['Number of reference samples (',num2str(nsamps),') should match with number of columns (',num2str(size(distance_cal,2)),') for distance_cal input'])
% 
%                 % Check number of beams dimensions match
%                 assert(isequal(size(Sv_cal,3),size(Echo_cal,3),size(distance_cal,3)),'Size of third dimension should match for Sv_cal, Echo_cal and distance_cal inputs')
%                 nbeams=size(Sv_cal,3);
% 
%                 % Check number of cells dimesnions match
%                 assert(isequal(size(Sv_cal,1),size(Echo_cal,1),size(distance_cal,1)),'Size of first dimension should match for Sv_cal, Echo_cal and distance_cal inputs')
% %                 ncells=size(Sv_cal,3);
%                 
%                 if nargout>1
%                     assert(isequal(size(Sv),size(Echo)), 'Sizes of Sv and Echo should match')
%                     assert(size(Sv,3)==nbeams, 'Size of third dimension of Sv should match with size of third dimension of Sv_cal');
%                 end
%                 % acoustics.SvSSC_Sassi.Kc_C get original or calibrated Kc and C values
%                 %
%                 %   Kc_C_out=Kc_C(obj,Sv_cal,Echo_cal,distance_cal) estimate the Kc and C
%                 %   given the volume backscatter strength, echo, distance from transducer.
%                 
%                 % Compute theoretical Sv
%                 ssc=[obj.reference_samples.mass_concentration]; % get sediment concentration of reference samples
%                 gsd=[obj.reference_samples.distribution];
%                 ks2=gsd.ks_squared(obj.transducer,[obj.reference_samples.sediment_density]); % get ks_squared for reference samples
%                 Sv_gsd=10*log10(ssc.*ks2); % Compute backscatter strength based on gsd
%                 
%                 % Find Sv measured near the samples
%                 cell_idx=acoustics.Sv2SSC_ConstantGSD.nearest_cells(distance_cal,obj.reference_samples);
%                 idx_4beams=sub2ind(size(Echo_cal),cell_idx,repmat(1:nsamps,[1,1,nbeams]),repmat(shiftdim(1:4,-1),[1, nsamps, 1]));
%                 Sv_cal_near=Sv_cal(idx_4beams);
%                 Echo_cal_near=Echo_cal(idx_4beams);
%                 
%                 %Do fit with Gostiaux type of equation, robust non linear fit
%                 Delta_sv=Sv_gsd-Sv_cal_near;
%                 Kc_C_out=nan(4,2);
%                 Sv_new=nan(size(Sv));
%                 for cb=1:4 % Fore each beam
%                     % Do the fitting
%                     Kc_C_out(cb,:)=nlinfit(...  % Make a non linear fit
%                         Echo_cal_near(1,:,cb),...     % x-values (independent)
%                         Delta_sv(1,:,cb),... % y-values (predicted)
%                         @(Kc_C,echo) acoustics.Sv2SSC_ConstantGSD.delta_Sv(Kc_C,obj.Kc_C_original, echo, obj.Er),...
%                         obj.Kc_C_original,... % Use original values as initial guess
%                         statset('RobustWgtFun','bisquare')); % use q robust bisquare weighting function
%                     % Compute new backscatter strength
%                     if nargout > 1
%                         Sv_new(:,cb)=... % new volume backscatter strength is:
%                             Sv(:,cb)+... % original + correction below
%                             acoustics.Sv2SSC_ConstantGSD.delta_Sv(Kc_C_out(cb,:), obj.Kc_C_original, Echo(:,cb), obj.Er);
%                     end
%                 end
%             else
%                 Kc_C_out=obj.Kc_C_original;
%                 if nargout>1
%                     Sv_new=Sv_cal;
%                 end
%             end
%         end
        function [b, gamma, stdb, stdgamma]=calibrate_b_gamma(obj)
            % obj
            sv=obj.sv_for_calibration();
            % calibration of b, for reference backscatter (paragraph 14)
            cell_idx=acoustics.Sv2SSC.nearest_cells(distance_ref,obj.reference_samples);
            Sv_ref=Sv_ref(cell_idx);
            beta_ref=acoustics.Sv2SSC.sv_to_beta(Sv_ref); % compute beta for ref (eqs. 12, Sassi et al.)
            k_ref=beta_ref./[obj.reference_samples.mass_concentration]; % compute k for ref (eqs. 12, Sassi et al.)
            if numel(k_ref)<3 % for less than two samples
                b=10*log10(k_ref)\Sv_ref; % simple linear fitting
            else % for more than two samples
                b=robustfit(10*log10(k_ref),Sv_ref,'bisquare',[],'off'); % least squares fit through origin (eq 14, Sassi et al.)... use robust fit here?
            end
            
            % calibration of specific attenuation (paragraph 15)
            cell_idx=acoustics.Sv2SSC.nearest_cells(distance_cal,obj.attenuation_samples);
            Sv_cal=Sv_cal(cell_idx);
            beta=acoustics.Sv2SSC.sv_to_beta(Sv_cal);
            ms=[obj.attenuation_samples.mass_concentration];
            k=beta_ref./ms;
            beta_all=acoustics.Sv2SSC_ConstantGSD.sv_to_beta(Sv_cal);
            int_beta=sum((beta_all(cell_idx(1):cell_idx(2)-1)+beta_all(cell_idx(1)+1:cell_idx(2))).*diff(distance_cal(cell_idx(1):cell_idx(2)))/2); % trapezoidal integration
            gamma=(k(1)-beta(2)./ms(2))./int_beta;
%             xi_e=5/log(10)*gamma;
        end
        function val=mass_concentration(Sv,b,gamma_e)
           %% sediment concentration calculation
            beta=acoustics.Sv2SSC_ConstantGSD.sv_to_beta(Sv);
            K_ref=beta(1,:,:).^b; % Kref from calibration

            int_beta=[zeros(1,size(beta,2), size(beta,3)); cumsum((beta(1:end-1,:,:)+beta(2:end,:,:)).*diff(d_sv,1,1)/2,1)];
            val=beta./(K_ref-gamma_e.*int_beta);
        end
    end
    methods(Static)
        function beta=sv_to_beta(sv)
            beta=10.^(sv/10);
        end
        function dsv=delta_Sv(Kc_C,Kc_C_original, Echo, Er)
            dsv=10*log10(...  % function to fit pars for
                (10.^(Kc_C(1).*Echo/10)-10.^(Kc_C(1).*Er/10))./...
                (10.^(Kc_C_original(1)*Echo/10)-10.^(Kc_C_original(1)*Er/10)))+...
                Kc_C(2)-Kc_C_original(2);
        end
    end
    
end


