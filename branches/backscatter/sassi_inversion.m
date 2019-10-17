function M=sassi_inversion(...
    sv_cal,...              backscatter strength
    d_sv_cal,...            depth of backscatter strength
    ref_samples,...     reference samples
    att_samples,...     attenuation samples
    sv,...
    d_sv...
    )%      time of attenuation samples
% Sediment concentration computation


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