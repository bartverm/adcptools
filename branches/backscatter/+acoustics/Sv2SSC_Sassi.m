classdef Sv2SSC_Sassi < handle
    properties
        Er(1,1) double {mustBeFinite, mustBeNonnegative} = 40
        reference_samples(:,1) acoustics.WaterSample
        attenuation_samples(:,1) acoustics.WaterSample
        piston_transducer(1,1) acoustics.PistonTransducer
        calibrate_density(1,1) logical = true
        calibrate_Kc_C(1,1) logical = true
        Kc_C_original(1,2) double = [0.45, -129.1]
        transducer(1,1) acoustics.PistonTransducer
    end
    properties(Dependent, SetAccess=private) % read only dependent properties
        all_samples
    end
    properties(Dependent)
        sediment_density(1,1) double
        Kc_C(1,2) double {mustBeFinite}
    end
    properties(Access=private)
        custom_sediment_density=2650;
    end
    methods % set and get methods
        function val=get.all_samples(obj)
            val=unique([obj.reference_samples(:); obj.attenuation_samples(:)]);
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
        function val=get.Kc_C(obj)
            if obj.calibrate_Kc_C
                ssc=[obj.reference_samples.mass_concentration];
                ks2=[obj.reference_samples.distribution.ks_squared(pt)];
                Sv=10*log10(ssc.*ks2);
                
                cell_idx=nearest_cells(d_sv_cal,ref_samples);
                sv_ref=sv_cal(cell_idx);

                % Using Deines type of equation
                Delta_sv=Sv-sv4(ind,:);
                x_maxi=ec4(ind,:);
                y_maxi=Delta_sv+kc*ec4(ind,:)+C;

                [pars, pars_maxi]=deal(nan(4,2));
                for cb=1:4
                    pars(cb,:)=robustfit(ec4(ind,cb),Delta_sv(:,cb));
                %     pars(cb,:)=[ones(numel(ind),1) ec4(ind,cb)]\Delta_sv(:,cb);
                %     pars_maxi(cb,:)=[ones(numel(ind),1), x_maxi(:,cb)]\y_maxi(:,cb);
                    pars_maxi(cb,:)=robustfit(x_maxi(:,cb),y_maxi(:,cb));
                end

                kc_maxi=pars_maxi(:,2); % Maxi adds .2 here... makes no sense to me
                C_maxi=pars_maxi(:,1);
                kc_this=kc+pars(:,2);
                C_this=C+pars(:,1);

                % Using Gostiaux type of equation
                pars_gost=nan(4,2);
                for cb=1:4
                    pars_gost(cb,:)=nlinfit(ec4(ind,cb),Delta_sv(:,cb),@(pars,x) 10*log10((10.^(pars(1).*x/10)-10.^(pars(1).*Er/10))./(10.^(kc*x/10)-10.^(kc*Er/10)))+pars(2)-C,[kc; C],statset('RobustWgtFun','bisquare'));
                end
                kc_gost=pars_gost(:,1);
                C_gost=pars_gost(:,2);
            else
                val=obj.Kc_C_original;
            end
        end
    end
    methods % ordinary methods

    end
    
    
end
    M=sassi_inversion(...
    sv_cal,...              backscatter strength
    d_sv_cal,...            depth of backscatter strength
    ref_samples,...         reference samples
    att_samples,...         attenuation samples
    sv,...
    d_sv,...
    varargin)
% Sediment concentration computation

%% Calibrate Kc and C (this is done with the ref_samples, near ADCP transducer)


%% calibration of b, for reference backscatter (paragraph 14)
cell_idx=nearest_cells(d_sv_cal,ref_samples);
sv_ref=sv_cal(cell_idx);
beta_ref=sv_to_beta(sv_ref); % compute beta for ref (eqs. 12, Sassi et al.)
k_ref=beta_ref./ref_samples.mass_concentration; % compute k for ref (eqs. 12, Sassi et al.)
if numel(k_ref)<3
    b=10*log10(k_ref)\sv_ref;
else
    b=robustfit(10*log10(k_ref),sv_ref,'bisquare',[],'off'); % least squares fit through origin (eq 14, Sassi et al.)... use robust fit here?
end
%% calibration of specific attenuation (paragraph 15)

% number of samples sets the number of attenuation values that can be
% obtained

cell_idx=nearest_cells(d_sv_cal,att_samples);
sv_att=sv_cal(cell_idx);
beta=sv_to_beta(sv_att);
ms=[att_samples(:).mass_concentration];
k=beta_ref./ms;
beta_all=sv_to_beta(sv_cal);
int_beta=sum((beta_all(cell_idx(1):cell_idx(2)-1)+beta_all(cell_idx(1)+1:cell_idx(2))).*diff(d_sv_cal(cell_idx(1):cell_idx(2)))/2); % trapezoidal integration

gamma_e=(k(1)-beta(2)./ms(2))./int_beta;
xi_e=5/log(10)*gamma_e;



%% sediment concentration calculation
beta=sv_to_beta(sv);
K_ref=beta(1,:,:).^b; % Kref from calibration

int_beta=[zeros(1,size(beta,2), size(beta,3)); cumsum((beta(1:end-1,:,:)+beta(2:end,:,:)).*diff(d_sv,1,1)/2,1)];
M=beta./(K_ref-gamma_e.*int_beta);

end


function beta=sv_to_beta(sv)
    beta=10.^(sv/10);
end
function idx=nearest_cells(dist, samples)
    d_sam=[samples(:).distance];
    [~,idx]=min((d_sam-dist).^2,[],1);
end