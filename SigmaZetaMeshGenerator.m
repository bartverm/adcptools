classdef SigmaZetaMeshGenerator < handle
% Base class to produce SigmaZetaMeshes
%
%   Subclasses should implement the get_mesh method.
%
%   SigmaZetaMeshGenerator methods:
%   get_mesh - returns the generated mesh
%
%   see also: Mesh, SigmaZetaMeshFromVMADCP
    methods(Abstract)
        mesh=get_mesh(obj)
    end
    
    methods(Static)
        function [nb, domain] = get_neighbors_and_domain(mesh)
            % provides neighbors of cells in mesh
            %
            % [nb, domain] = get_neighbors_and_domain(mesh) returns the 
            % % array nb with size 4 x ncells holding the index of the 
            % % neighbor to right, top left and bottom in the first to 
            % fourth row respectively.
            %
            %                     | 2 -> top    |
            %           __________|_____________|___________
            %                     |             |
            %           3 -> left |    cell     | 1 -> right
            %           __________|_____________|___________
            %                     |             |
            %                     | 4 -> bottom | 
            %
            % The domain array hold an index for each cell depending on
            % where the cell is located in the mesh:
            %
            %               4 |         3         | 2
            %              ___|___________________|___ 
            %                 |                   |
            %               5 |         0         | 1
            %              ___|___________________|___
            %                 |                   |
            %               6 |         7         | 8
            %
            %   9: degenerate cell, i.e. having outer border on more than
            %   two sides, or having it on the bottom and top or left and
            %   right
            %

            nb = nan(mesh.ncells,4);
            sig = reshape(mesh.sig_center,1,[]);
            Sig = nan(size(mesh.col_to_mat));
            cell_to_mat = mesh.cell_to_mat;
            idx_mat = nan(size(Sig));
            idx_mat(cell_to_mat) = 1:mesh.ncells;
            min_idx = min(idx_mat,[],1,'omitnan');
            max_idx = max(idx_mat,[],1,'omitnan');
            Sig(cell_to_mat) = sig;
            isnan_mat = isnan(Sig);
            dsig = mean(diff(Sig,1,1),1,'omitnan');
            Sig = cumsum([Sig(1,:); repmat(dsig,[size(Sig,1)-1 1])],1);
            idx_mat = cumsum([idx_mat(1,:); ones(size(Sig)-[1 0])]);

            % find right neighbor
            dSig_right = diff(Sig,1,2);
            dsig = abs(dsig);
            adjustIdx_right = floor(dSig_right./dsig(2:end)+.5); 
            idx_right = idx_mat(:,2:end) + adjustIdx_right;
            idx_right = min(idx_right, max_idx(2:end));
            idx_right = max(idx_right, min_idx(2:end));
            idx_right = [idx_right nan(size(idx_right,1),1)];
            idx_right(isnan_mat) = nan;
            nb(:,1) = idx_right(~isnan_mat);

            % find top neighbor
            idx_top = idx_mat(1:end-1,:);
            idx_top = [nan(1,size(idx_top,2)); idx_top];
            idx_top(isnan_mat) = nan;
            nb(:,2) = idx_top(~isnan_mat);

            % find left neighbor
            dSig_left = -dSig_right;
            adjustIdx_left = floor(dSig_left./dsig(1:end-1)+.5); 
            idx_left = idx_mat(:,1:end-1) + adjustIdx_left;
            idx_left = min(idx_left, max_idx(1:end-1));
            idx_left = max(idx_left, min_idx(1:end-1));
            idx_left = [nan(size(idx_left,1),1) idx_left];
            idx_left(isnan_mat) = nan;
            nb(:,3) = idx_left(~isnan_mat);

            % find bottom neighbor
            idx_bot = idx_mat;
            idx_bot(isnan_mat) = nan;
            idx_bot = [idx_bot; nan(1,size(idx_bot,2))];
            idx_bot = idx_bot(2:end,:);
            nb(:,4) = idx_bot(~isnan_mat);

            % find domains if needed
            if nargout < 2
                return
            end
            domain = ones(mesh.ncells,1)*9;
            nbnan = isnan(nb);
            nnan = sum(nbnan,2);
            is_corn = nnan == 2;
            is_side = nnan ==1;
            domain(nnan == 0) = 0;
            domain(is_side & nbnan(:,1)) = 1;
            domain(is_side & nbnan(:,2)) = 3;
            domain(is_side & nbnan(:,3)) = 5;
            domain(is_side & nbnan(:,4)) = 7;
            domain(is_corn & all(nbnan(:,[1 2]),2)) = 2;
            domain(is_corn & all(nbnan(:,[2 3]),2)) = 4;
            domain(is_corn & all(nbnan(:,[3 4]),2)) = 6;
            domain(is_corn & all(nbnan(:,[4 1]),2)) = 8;
        end
        function jac = get_jacobian(mesh)
            jac = nan(mesh.ncells,2,2);
            n = reshape(mesh.n_middle(mesh.col_to_cell),[],1);
            sig = mesh.sig_center;
            has_left_right = ismember(mesh.domains, [0 3 7]);
            has_top_bottom = ismember(mesh.domains, [0 1 5]);
            nb = mesh.neighbors;
            jac(has_left_right, 1, 1) = .5 * (...
                n(nb(has_left_right, 3))-...
                n(nb(has_left_right, 1))); %dn/di
            jac(has_top_bottom, 1, 2) = 0; %dn/dj
            jac(has_left_right, 2, 1) = .5 * (...
                sig(nb(has_left_right, 3))-...
                sig(nb(has_left_right, 1))); % dsig/di
            jac(has_top_bottom, 2, 2) = .5 * (...
                sig(nb(has_top_bottom, 2))-...
                sig(nb(has_top_bottom, 4))); % dsig/dj
        end
    end
end