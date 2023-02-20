function C1 = assembleC1(LBS)

% Function that assembles cell-based continuity equation
cmesh = LBS.mesh;
ZB = LBS.bathy;
T = LBS.velocity_model;
eta = LBS.adcp.water_level_object;
B = 1;
L = 1;
dx = 2; dy = 2;
t = LBS.mesh.time;

continuity_possible = 0;
if (T.s_order(1) > 0) && (T.n_order(2) > 0) && (T.sigma_order(3) > 0 || T.z_order(3) > 0)
    continuity_possible = 1;
end



for idx = 1:cmesh.ncells
sig = cmesh.sig_center(idx);
D0 = eta.parameters(1) - LBS.mesh.zb_middle(LBS.mesh.col_to_cell(idx));
D0x = -1/(2*dx)*(ZB.get_depth(cmesh.x_middle(LBS.mesh.col_to_cell(idx))+dx,cmesh.y_middle(LBS.mesh.col_to_cell(idx)),t)-...
    ZB.get_depth(cmesh.x_middle(LBS.mesh.col_to_cell(idx))-dx,cmesh.y_middle(LBS.mesh.col_to_cell(idx)),t));% Assumption of zero gradient in x-direction
D0y = -1/(2*dy)*(ZB.get_depth(cmesh.x_middle(LBS.mesh.col_to_cell(idx)),cmesh.y_middle(LBS.mesh.col_to_cell(idx))+dy,t)-...
    ZB.get_depth(cmesh.x_middle(LBS.mesh.col_to_cell(idx)),cmesh.y_middle(LBS.mesh.col_to_cell(idx))-dy,t));% Assumption of zero gradient in x-direction

Cj{idx} = zeros([1+2*numel(eta.constituents), sum(T.npars)]);
% ind = zeros([3, numel(names)]);
% First: subtidal equation

if continuity_possible
    ind0 = [find(strcmp(T.names, 'd^1u/dx^1: M0A')) , ...
        find(strcmp(T.names, 'd^1u/dsig^1: M0A')) , ...
        find(strcmp(T.names, 'd^1v/dy^1: M0A')) , ...
        find(strcmp(T.names, 'd^1v/dsig^1: M0A')) , ...
        find(strcmp(T.names, 'd^1w/dsig^1: M0A'))];

    terms0 = [B*D0, B*(1-sig)*D0x, L*D0, L*(1-sig)*D0y, L*B];

    Cj{idx}(1,ind0) = terms0;


    % Second: Tidal equations for each constituent
    for i=1:numel(T.constituentsU) %Fundamental choice: All constituents are the same!!
        const = T.constituentsU{i};

        ind{i,1} = [find(strcmp(T.names, sprintf('%s%s%s', 'd^1u/dx^1: ', const,'A'))) , ...
            find(strcmp(T.names, 'd^1u/dx^1: M0A')) , ...
            find(strcmp(T.names, sprintf('%s%s%s', 'd^1u/dsig^1: ', const,'A'))) ,...
            find(strcmp(T.names, sprintf('%s%s%s', 'd^1v/dy^1: ', const,'A'))) , ...
            find(strcmp(T.names, 'd^1v/dy^1: M0A')) , ...
            find(strcmp(T.names, sprintf('%s%s%s', 'd^1v/dsig^1: ', const,'A'))) , ...
            find(strcmp(T.names, sprintf('%s%s%s', 'd^1w/dsig^1: ', const,'A')))];

        terms{i,1} = [B*D0, B*eta.parameters(2*i), B*(1-sig)*D0x, L*D0, B*eta.parameters(2*i), L*(1-sig)*D0y, L*B];
        Cj{idx}(2*i,ind{i,1}) = terms{i,1};


        ind{i,2} = [find(strcmp(T.names, sprintf('%s%s%s', 'd^1u/dx^1: ', const,'B'))) , ...
            find(strcmp(T.names, 'd^1u/dx^1: M0A')) , ...
            find(strcmp(T.names, sprintf('%s%s%s', 'd^1u/dsig^1: ', const,'B'))) , ...
            find(strcmp(T.names, sprintf('%s%s%s', 'd^1v/dy^1: ', const,'B'))) , ...
            find(strcmp(T.names, 'd^1v/dy^1: M0A')) , ...
            find(strcmp(T.names, sprintf('%s%s%s', 'd^1v/dsig^1: ', const,'B'))) , ...
            find(strcmp(T.names, sprintf('%s%s%s', 'd^1w/dsig^1: ', const,'B')))];

        terms{i,2} = [B*D0, B*eta.parameters(2*i+1), B*(1-sig)*D0x, L*D0, B*eta.parameters(2*i+1), L*(1-sig)*D0y, L*B];
        Cj{idx}(2*i+1,ind{i,2}) = terms{i,2};
    end
else
    disp('Warning: No continuity matrix assembled')
end

end
C1 = spblkdiag(Cj{:});
% To make sure it has the right dimensions
% C(cmesh.ncells*(1+2*numel(eta.constituents)), cmesh.ncells*sum(T.npars)) = 0;