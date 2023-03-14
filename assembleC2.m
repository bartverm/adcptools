function Cg = assembleC2(LBS)
par_names_tot = repmat(LBS.velocity_model.names,1,LBS.mesh.ncells);
Np = sum(LBS.velocity_model.npars);
idx = 1;
for j = 1:LBS.mesh.ncells % loop trough every cell
    for i=1:Np
        par_names_tot{1,idx} = sprintf('cell %i: %s',j, par_names_tot{1,idx});
        idx = idx + 1;
    end
    [neighbors(j,:), dom(j)] = LBS.mesh.get_neighbors(j);
end
row_idx = 1;
eta = LBS.adcp.water_level_object;
Nt = length(eta.names);
D0x = 0;
D0y = 0;
sig_center = LBS.mesh.sig_center;
n_center = LBS.mesh.n_middle(LBS.mesh.col_to_cell);
Bt = 1;
Lt = 1;
rows = []; cols = []; vals = [];
for c = 1:LBS.mesh.ncells %rows
    adj = neighbors(c,:);
    doma = dom(c);
    D0 = LBS.adcp.water_level_object.parameters(1) - LBS.mesh.zb_middle(LBS.mesh.col_to_cell(c));

    if doma==0 || doma==1 || doma == 5
        dsig = sig_center(adj(2))-sig_center(adj(4)); % central difference
    elseif doma ==2 || doma==3|| doma==4
        dsig = sig_center(c)-sig_center(adj(4)); % one-sided difference
    else
        dsig = sig_center(adj(2))-sig_center(c); % one-sided difference
    end

    if doma==0 || doma==3 || doma == 7
        dn = n_center(adj(1))-n_center(adj(3)); % central difference
    elseif doma ==1 || doma==2|| doma==8
        dn = n_center(c)-n_center(adj(3)); % one-sided difference
    else
        dn = n_center(adj(1))-n_center(c); % one-sided difference
    end
    % Subtidal
    row = row_idx*ones(1,9); % 9 terms
    val = [Bt*D0, Bt*(1-sig_center(c))*D0x/dsig, -Bt*(1-sig_center(c))*D0x/dsig,...
        Lt*D0/dn, -Lt*D0/dn, Lt*(1-sig_center(c))*D0y/dsig, -Lt*(1-sig_center(c))*D0y/dsig,...
        Lt*Bt/dsig, -Lt*Bt/dsig];
    [row, col, val] = dom2rowcol(par_names_tot, doma, c, adj, 'M0A', row, val);
    rows = [rows row];
    cols = [cols col];
    vals = [vals val];
    row_idx = row_idx + 1;
    if length(col)==length(val)
        disp('')
    else
        disp('er')
    end
    for co = 2:Nt
        nam = LBS.adcp.water_level_object.names{co}(end-2:end);
        row = row_idx*ones(1,11);
        val = [Bt*D0, Bt*(1-sig_center(c))*D0x/dsig, -Bt*(1-sig_center(c))*D0x/dsig,...
             Lt*D0/dn, -Lt*D0/dn, Lt*(1-sig_center(c))*D0y/dsig, -Lt*(1-sig_center(c))*D0y/dsig,...
             Lt*Bt/dsig, -Lt*Bt/dsig, Bt*eta.parameters(co), Lt*eta.parameters(co)];
        [row, col, val] = dom2rowcol(par_names_tot, doma, c, adj, nam, row, val);
        rows = [rows row];
        cols = [cols col];
        vals = [vals val];
        row_idx = row_idx + 1;
    end
end
Cg = sparse(rows, cols, vals, LBS.mesh.ncells*(2*length(LBS.velocity_model.constituentsU)+1), LBS.mesh.ncells*sum(LBS.velocity_model.npars));
% spy(Cg)
end