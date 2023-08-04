function [tid_pars, tid_names] = ab2Aphi(pars, names)
% Function that converts tidal constituents and gradients to A - phi
% representation
tid_pars = zeros(size(pars));
tid_names = names;
ind = 1;
while ind <= length(names)
    if contains(names{ind}, 'M0') % Subtidal, spatially constant and varying
        tid_pars(:,ind) = pars(:,ind);
        tid_names{1, ind} = strrep(names{1,ind}, 'M0A', 'M0');
        ind = ind + 1;
    elseif ~contains(names{ind}, 'd') % Tidal, spatially constant
        tid_pars(:,ind) = sqrt(pars(:,ind).^2 + pars(:,ind+1).^2);
        tid_pars(:,ind+1) = atan2(pars(:,ind+1),pars(:,ind));
        tid_names{1, ind} = strrep(names{1,ind}, 'A', ': A');
        tid_names{1, ind+1} = strrep(names{1,ind+1}, 'B', ': \phi');
        ind = ind + 2;
        tempa = pars(:,ind); tempb = pars(:,ind+1);
    else % Tidal, spatially varying -> Chain rule
        tid_pars(:,ind) = (tempa.*pars(:,ind) + tempb.*pars(:,ind+1))./(sqrt(tempa.^2 + tempb.^2));
        tid_pars(:,ind + 1) = (tempa.*pars(:,ind+1) - tempb.*pars(:,ind))./(tempa.^2 + tempb.^2);
        tid_names{1, ind} = strrep(names{1,ind}, 'A', ': A');
        tid_names{1, ind+1} = strrep(names{1,ind+1}, 'B', ': \phi');        
        ind = ind+2;
    end

end