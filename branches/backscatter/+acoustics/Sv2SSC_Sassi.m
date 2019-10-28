classdef Sv2SSC_Sassi < handle
% Sv2SSC_Sassi computes mass conc. from backscatter (Sassi et al., 2012)
%
%   obj=Sv2SSC_Sassi() Constructs default object
%
%   Sv2SSC_Sassi properties:
%   Er - background noise
%   reference_samples - water samples for reference, near transducer
%   attenuation_samples - water samples for atten. calibr in water column
%   calibrate_density - whether density of particles should be estimated
%   sediment_density - calibrated or fixed sediment density
%   calibrate_Kc_C - whether to calibrate the Kc and C instrument constants
%   Kc_C_original - initial values of Kc and C
%   transducer - piston transducer object
%
%   Sv2SSC_Sassi read only property:
%   all_samples - returns all unique water samples
%
%   Sv2SSC_Sassi methods:
%   Kc_C - returns the original or calibrated Kc and C values
%   b_gamma - calibrate the b and gamma calibration constants
%   mass_concentration - estimate the mass concentration profiles
%
%   Sv2SSC_Sassi static methods:
%   sv_to_beta - compute beta from Sv values
%   nearest_cells - find cells nearest to sample depth
%   delta_Sv - find correction for Sv based on new Kc and C
%
%   see also: svADCP

% TODO: work with ADCP objects, make Kc_C a dependent property similar to
% sediment_density, testing!
    properties
        % acoustics.Sv2SSC_Sassi/Er property
        %
        % scalar, non-negative, finite double indicating background noise
        % echo in counts. Default: 40.
        %
        % see also: acoustics.Sv2SSC_Sassi
        Er(1,1) double {mustBeFinite, mustBeNonnegative} = 40
        
        % acoustics.Sv2SSC_Sassi/reference_samples property
        %
        % Vector of acoustic.WaterSample objects representing the reference
        % samples collected near the instruments' transducer.
        %
        % see also: acoustics.Sv2SSC_Sassi
        reference_samples(:,1) acoustics.WaterSample
        
        % acoustics.Sv2SSC_Sassi/attenuation_samples property
        %
        % Vector of acoustic.WaterSample objects representing the samples 
        % used for calibration of the attenuation (gamma) and 
        % backscatter (b) calibration constants
        %
        % see also: acoustics.Sv2SSC_Sassi
        attenuation_samples(:,1) acoustics.WaterSample
        
        % acoustics.Sv2SSC_Sassi/calibrate_density property
        %
        % scalar, logical value indicating whether the density of the
        % particles should be calibrated from the water samples
        %
        % see also: acoustics.Sv2SSC_Sassi
        calibrate_density(1,1) logical = true
        
        % acoustics.Sv2SSC_Sassi/calibrate_Kc_C property
        %
        % scalar, logical value indicating wheter the Kc and C instrument
        % constants should be estimated based on grain size distribution
        % and concentrations of reference water samples
        %
        % see also: acoustics.Sv2SSC_Sassi
        calibrate_Kc_C(1,1) logical = true
        
        % acoustics.Sv2SSC_Sassi/Kc_C_original property
        %
        % two element vector with double values holding the Kc and C
        % constants originally used for the volume backscatter computations
        %
        % see also: acoustics.Sv2SSC_Sassi, svADCP
        Kc_C_original(1,2) double = [0.45, -129.1]
        
        % acoustics.Sv2SSC_Sassi/transducer property
        %
        % scalar acoustic.PistonTransducer object needed to estimate sound
        % properties
        %
        % see also: acoustics.Sv2SSC_Sassi
        transducer(1,1) acoustics.PistonTransducer
    end
    properties(Dependent, SetAccess=private)
        % acoustics.Sv2SSC_Sassi/all_samples read only property
        %
        % holds all unique water samples taken from the reference and
        % attenuation samples
        %
        % see also: acoustics.Sv2SSC_Sassi
        all_samples
    end
    properties(Dependent)
        % acoustics.Sv2SSC_Sassi/sediment_density property
        %
        % a custom sediment density can be assigned to this property, which
        % by default is 2650. If calibration of sediment density is enabled
        % the value will return the calibrated sediment density. Assigning
        % to this variable will disable the sediment density calibration
        % and a warning is generated.
        %
        % see also: acoustics.Sv2SSC_Sassi
        sediment_density(1,1) double
    end
    properties(Access=private)
        custom_sediment_density=2650;
    end
    methods 
        %%% Set and get methods %%%
        function val=get.all_samples(obj)
            val=unique([obj.reference_samples; obj.attenuation_samples]);
        end
        function val=get.sediment_density(obj)
            if obj.calibrate_density
                val=[obj.all_samples.sediment_density];
                val(~isfinite(val))=[];
                val=mean(val);
            else
                val=obj.custom_sediment_density;              
            end
        end
        function set.sediment_density(obj,val)
            if obj.calibrate_density
                warning('Turning off density calibration')
                obj.calibrate_density=false;
            end
            obj.custom_sediment_density=val;
        end
        
        %%% Ordinary methods %%%
        function [Kc_C_out, Sv_new]=Kc_C(obj,Sv_cal,Echo_cal,distance_cal,Sv, Echo)
            if obj.calibrate_Kc_C % if we want the calibration

                % check number of samples dimensions match
                samps=obj.reference_samples;
                nsamps=numel(samps);
                assert(size(Sv_cal,2)==nsamps, ['Number of reference samples (',num2str(nsamps),') should match with number of columns (',num2str(size(Sv_cal,2)),') for Sv_cal input'])
                assert(size(Echo_cal,2)==nsamps, ['Number of reference samples (',num2str(nsamps),') should match with number of columns (',num2str(size(Echo_cal,2)),') for Echo_cal input'])
                assert(size(distance_cal,2)==nsamps, ['Number of reference samples (',num2str(nsamps),') should match with number of columns (',num2str(size(distance_cal,2)),') for distance_cal input'])

                % Check number of beams dimensions match
                assert(isequal(size(Sv_cal,3),size(Echo_cal,3),size(distance_cal,3)),'Size of third dimension should match for Sv_cal, Echo_cal and distance_cal inputs')
                nbeams=size(Sv_cal,3);

                % Check number of cells dimesnions match
                assert(isequal(size(Sv_cal,1),size(Echo_cal,1),size(distance_cal,1)),'Size of first dimension should match for Sv_cal, Echo_cal and distance_cal inputs')
%                 ncells=size(Sv_cal,3);
                
                if nargout>1
                    assert(isequal(size(Sv),size(Echo)), 'Sizes of Sv and Echo should match')
                    assert(size(Sv,3)==nbeams, 'Size of third dimension of Sv should match with size of third dimension of Sv_cal');
                end
                % acoustics.SvSSC_Sassi.Kc_C get original or calibrated Kc and C values
                %
                %   Kc_C_out=Kc_C(obj,Sv_cal,Echo_cal,distance_cal) estimate the Kc and C
                %   given the volume backscatter strength, echo, distance from transducer.
                
                % Compute theoretical Sv
                ssc=[obj.reference_samples.mass_concentration]; % get sediment concentration of reference samples
                gsd=[obj.reference_samples.distribution];
                ks2=gsd.ks_squared(obj.transducer,[obj.reference_samples.sediment_density]); % get ks_squared for reference samples
                Sv_gsd=10*log10(ssc.*ks2); % Compute backscatter strength based on gsd
                
                % Find Sv measured near the samples
                cell_idx=acoustics.Sv2SSC_Sassi.nearest_cells(distance_cal,obj.reference_samples);
                idx_4beams=sub2ind(size(Echo_cal),cell_idx,repmat(1:nsamps,[1,1,nbeams]),repmat(shiftdim(1:4,-1),[1, nsamps, 1]));
                Sv_cal_near=Sv_cal(idx_4beams);
                Echo_cal_near=Echo_cal(idx_4beams);
                
                %Do fit with Gostiaux type of equation, robust non linear fit
                Delta_sv=Sv_gsd-Sv_cal_near;
                Kc_C_out=nan(4,2);
                Sv_new=nan(size(Sv));
                for cb=1:4 % Fore each beam
                    % Do the fitting
                    Kc_C_out(cb,:)=nlinfit(...  % Make a non linear fit
                        Echo_cal_near(1,:,cb),...     % x-values (independent)
                        Delta_sv(1,:,cb),... % y-values (predicted)
                        @(Kc_C,echo) acoustics.Sv2SSC_Sassi.delta_Sv(Kc_C,obj.Kc_C_original, echo, obj.Er),...
                        obj.Kc_C_original,... % Use original values as initial guess
                        statset('RobustWgtFun','bisquare')); % use q robust bisquare weighting function
                    % Compute new backscatter strength
                    if nargout > 1
                        Sv_new(:,cb)=... % new volume backscatter strength is:
                            Sv(:,cb)+... % original + correction below
                            acoustics.Sv2SSC_Sassi.delta_Sv(Kc_C_out(cb,:), obj.Kc_C_original, Echo(:,cb), obj.Er);
                    end
                end
            else
                Kc_C_out=obj.Kc_C_original;
                if nargout>1
                    Sv_new=Sv_cal;
                end
            end
        end
        function [b, gamma]=b_gamma(obj, Sv_ref, distance_ref, Sv_cal, distance_cal)
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
            beta_all=acoustics.Sv2SSC_Sassi.sv_to_beta(Sv_cal);
            int_beta=sum((beta_all(cell_idx(1):cell_idx(2)-1)+beta_all(cell_idx(1)+1:cell_idx(2))).*diff(distance_cal(cell_idx(1):cell_idx(2)))/2); % trapezoidal integration
            gamma=(k(1)-beta(2)./ms(2))./int_beta;
%             xi_e=5/log(10)*gamma;
        end
        function val=mass_concentration(Sv,b,gamma_e)
           %% sediment concentration calculation
            beta=acoustics.Sv2SSC_Sassi.sv_to_beta(Sv);
            K_ref=beta(1,:,:).^b; % Kref from calibration

            int_beta=[zeros(1,size(beta,2), size(beta,3)); cumsum((beta(1:end-1,:,:)+beta(2:end,:,:)).*diff(d_sv,1,1)/2,1)];
            val=beta./(K_ref-gamma_e.*int_beta);
        end
    end
    methods(Static)
        function beta=sv_to_beta(sv)
            beta=10.^(sv/10);
        end
        function idx=nearest_cells(dist, samples)
            d_sam=[samples(:).distance];
            [~,idx]=min((d_sam-dist).^2,[],1);
        end
        function dsv=delta_Sv(Kc_C,Kc_C_original, Echo, Er)
            dsv=10*log10(...  % function to fit pars for
                (10.^(Kc_C(1).*Echo/10)-10.^(Kc_C(1).*Er/10))./...
                (10.^(Kc_C_original(1)*Echo/10)-10.^(Kc_C_original(1)*Er/10)))+...
                Kc_C(2)-Kc_C_original(2);
        end
    end
    
end


