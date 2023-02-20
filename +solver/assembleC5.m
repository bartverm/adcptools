function [C, rhsvec] = assembleC5(LBS)

% Function that assembles cell-based kinematic boundary conditions.

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

cmesh = LBS.mesh;
ZB = LBS.bathy;
T = LBS.velocity_model;
eta = LBS.adcp.water_level_object;
eta.get_omega();
B = 1;
L = 1;
dx = 2; dy = 2;
t = LBS.mesh.time;

% Compute relevant indices

ind0 = [find(strcmp(T.names, sprintf('%s%s%s', 'u0: M0A'))) , ...
    find(strcmp(T.names, sprintf('%s%s%s', 'd^1u/dsig^1: M0A'))) ,...
    find(strcmp(T.names, sprintf('%s%s%s', 'v0: M0A'))) , ...
    find(strcmp(T.names, sprintf('%s%s%s', 'd^1v/dsig^1: M0A'))) , ...
    find(strcmp(T.names, sprintf('%s%s%s', 'w0: M0A'))) , ...
    find(strcmp(T.names, sprintf('%s%s%s', 'd^1w/dsig^1: M0A')))];

for i=1:numel(T.constituentsU) %Fundamental choice: All constituents are the same!!
    const = T.constituentsU{i};
    ind{i,1} = [find(strcmp(T.names, sprintf('%s%s%s', 'u0: ', const,'A'))) , ...
        find(strcmp(T.names, sprintf('%s%s%s', 'd^1u/dsig^1: ', const,'A'))) ,...
        find(strcmp(T.names, sprintf('%s%s%s', 'v0: ', const,'A'))) , ...
        find(strcmp(T.names, sprintf('%s%s%s', 'd^1v/dsig^1: ', const,'A'))) , ...
        find(strcmp(T.names, sprintf('%s%s%s', 'w0: ', const,'A'))) , ...
        find(strcmp(T.names, sprintf('%s%s%s', 'd^1w/dsig^1: ', const,'A')))];

    ind{i,2} = [find(strcmp(T.names, sprintf('%s%s%s', 'u0: ', const,'B'))) , ...
        find(strcmp(T.names, sprintf('%s%s%s', 'd^1u/dsig^1: ', const,'B'))) ,...
        find(strcmp(T.names, sprintf('%s%s%s', 'v0: ', const,'B'))) , ...
        find(strcmp(T.names, sprintf('%s%s%s', 'd^1v/dsig^1: ', const,'B'))) , ...
        find(strcmp(T.names, sprintf('%s%s%s', 'w0: ', const,'B'))) , ...
        find(strcmp(T.names, sprintf('%s%s%s', 'd^1w/dsig^1: ', const,'B')))];
end

% Now, the indices of the relevant terms are known. Proceed to fill in
el_mat_size  = [1+2*numel(eta.constituents), sum(T.npars)];
for idx = 1:cmesh.ncells
    adj = neighbors(idx,:);
    doma = dom(idx);
    Cj{idx} = zeros(el_mat_size);
    rhs{idx} = zeros([el_mat_size(1),1]);
    if doma==0 || doma==1 || doma == 5 % Internal cells (in terms of sigma)
        % Do nothing
    else
        bot = (doma ==6 || doma==7|| doma==8);
        surf = (doma ==2 || doma==3|| doma==4); %Surface cells
        if surf
            dsig = 1-cmesh.sig_center(idx); % per definition of sigma. Positive
            % ind = zeros([3, numel(names)]);
            % First: subtidal equation
            terms0 = [0, 0, 0, 0, 1, dsig];
            Cj{idx}(1,ind0) = terms0;
            rhs{idx}(1,1) = 0;
            for i=1:numel(T.constituentsU) %Fundamental choice: All constituents are the same!!
                terms{i,1} = terms0;
                Cj{idx}(2*i,ind{i,1}) = terms{i,1};
                terms{i,2} = terms0;
                Cj{idx}(2*i+1,ind{i,2}) = terms{i,2};
                rhs{idx}(2*i, 1) = eta.omega(i)*eta.parameters(2*i);
                rhs{idx}(2*i+1, 1) = -eta.omega(i)*eta.parameters(2*i+1);
            end
        elseif bot % Bottom cells
            dsig = cmesh.sig_center(idx); % 
            % Assemble element matrix
%             sig = cmesh.sig_center(idx);
            %             D0 = eta.parameters(1) - LBS.mesh.zb_middle(LBS.mesh.col_to_cell(idx));
            zbx = 1/(2*dx)*(ZB.get_depth(cmesh.x_middle(LBS.mesh.col_to_cell(idx))+dx,cmesh.y_middle(LBS.mesh.col_to_cell(idx)),t)-...
                ZB.get_depth(cmesh.x_middle(LBS.mesh.col_to_cell(idx))-dx,cmesh.y_middle(LBS.mesh.col_to_cell(idx)),t));% Assumption of zero gradient in x-direction
            zby = 1/(2*dy)*(ZB.get_depth(cmesh.x_middle(LBS.mesh.col_to_cell(idx)),cmesh.y_middle(LBS.mesh.col_to_cell(idx))+dy,t)-...
                ZB.get_depth(cmesh.x_middle(LBS.mesh.col_to_cell(idx)),cmesh.y_middle(LBS.mesh.col_to_cell(idx))-dy,t));% Assumption of zero gradient in x-direction

            terms0 = [zbx, -dsig*zbx, zby, -dsig*zby, -1, dsig];
            Cj{idx}(1,ind0) = terms0;
            for i=1:numel(T.constituentsU) %Fundamental choice: All constituents are the same!!
                terms{i,1} = terms0;
                Cj{idx}(2*i,ind{i,1}) = terms{i,1};
                terms{i,2} = terms0;
                Cj{idx}(2*i+1,ind{i,2}) = terms{i,2};
            end
        end
    end
end


C = spblkdiag(Cj{:});
rhsvec = cell2mat(rhs');
% To make sure it has the right dimensions
% C(cmesh.ncells*(1+2*numel(eta.constituents)), cmesh.ncells*sum(T.npars)) = 0;