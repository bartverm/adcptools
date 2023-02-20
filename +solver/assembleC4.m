function C4 = assembleC4(LBS)
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
sig_center = LBS.mesh.sig_center;
n_center = LBS.mesh.n_middle(LBS.mesh.col_to_cell);
rows = []; cols = []; vals = [];
for c = 1:LBS.mesh.ncells %rows
    adj = neighbors(c,:);
    doma = dom(c);
    if doma==0 || doma==1 || doma == 5 % Internal cells (in terms of sigma)
        dsig = sig_center(adj(2))-sig_center(adj(4)); % central difference
    elseif doma ==2 || doma==3|| doma==4 %Surface cells
        dsig = sig_center(c)-sig_center(adj(4)); % one-sided difference
    else % Bottom cells
        dsig = sig_center(adj(2))-sig_center(c); % one-sided difference
    end

    if doma==0 || doma==3 || doma == 7 % Internal cells (laterally)
        dn = n_center(adj(1))-n_center(adj(3)); % central difference
    elseif doma ==1 || doma==2|| doma==8 % Right side of domain
        dn = n_center(c)-n_center(adj(3)); % one-sided difference
    else % Left side of domain
        dn = n_center(adj(1))-n_center(c); % one-sided difference
    end
    for co = 1:Nt
        nam = LBS.adcp.water_level_object.names{co}(end-2:end);
        [row, col, val] = dom2rowcol2(par_names_tot, doma, c, adj, nam, dn, dsig);
        rows = [rows row_idx+row];
        cols = [cols col];
        vals = [vals val];
        row_idx = max(rows) + 1;
    end
end
C4 = sparse(rows, cols, vals, 6*LBS.mesh.ncells*(2*length(LBS.velocity_model.constituentsU)+1), LBS.mesh.ncells*sum(LBS.velocity_model.npars));

end