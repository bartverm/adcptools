classdef SigmaZetaMesh < Mesh
% Defines a sigma-z mesh
%   
%   The SigmaZetaMesh should be generated with a SigmaZetaMeshGenerator.
%
%   The mesh consists of verticals. Each vertical has a certain number of
%   cells. In each vertical the cells follow the bed. The name Sigma-Zeta
%   refers to the fact that this mesh is a hyrbid between a z-mesh and a
%   sigma-mesh. It tries to combine the best of both worlds: it follows
%   nicely the bed just like a sigma mesh would, but the vertical spacing
%   of the cells is constant, just like a z-mesh would. This means that
%   each cell will hold approximately the same number of adcp velocity
%   estimates.
%
%   Each mesh cell consists of six edges, called:
%           *         ---> top-middle
%         /   \
%       /       \
%     /          *    ---> top-right
%   *             |   ---> top-left 
%   |       *     |   ---> bottom-middle
%   |     /   \   |
%   |   /       \ |
%   | /           *   ---> bottom-right
%   *                 ---> bottom-left
%
%   Note that the left, middle and right edges are always vertically
%   stacked and share the same n-coordinate.
%
%   Data can be stored either in a vector with a value for each cell or in
%   a matrix which has a toplogy similar to the mesh. The indexing
%   properties (see below) help map between these formats.
%
%   SigmaZetaMesh properties (read only):
%   * General
%   xs - defines the cross-section of the mesh
%   water_level - the water level for the mesh
%   ncells - number of cells in the mesh
%   
%   * Edge Position
%   z_bottom_left - z-coordinate of the bottom-left edge
%   z_top_left - z-coordinate of the top-left edge
%   z_bottom_mid - z-coordinate of the bottom-middle edge
%   z_top_mid - z-coordinate of the top-middle edge
%   z_bottom_right - z-coordinate of the bottom-right edge
%   z_top_right - z-coordinate of the top-right edge
%   sig_bottom_left - sigma-coordinate of the bottom-left edge
%   sig_top_left - sigma-coordinate of the top-left edge
%   sig_bottom_mid - sigma-coordinate of the bottom-middle edge
%   sig_top_mid - sigma-coordinate of the top-middle edge
%   sig_bottom_right - sigma-coordinate of the bottom-right edge
%   sig_top_right - sigma-coordinate of the top-right edge
%   n_left - n-coordinate of the left edges
%   n_middle - n-coordinate of the middle edges
%   n_right - n-coordinate of the right edges
%   x_left - x-coordinate of the left edges
%   x_middle - x-coordinate of the middle edges
%   x_right - x-coordinate of the right edges
%   y_left - y-coordinate of the left edges
%   y_middle - y-coordinate of the middle edges
%   y_right - y-coordinate of the right edges
%
%   * Bed position
%   zb_left - z-coordinate of bed at left edges
%   zb_middle - z-coordinate of bed at middle edges
%   zb_right - z-coordinate of bed at right edges
%   zb_all - z-coordinate of bed at all edges from left to right
%
%   * Water surface position
%   nw - n-coordinates of water surface boundaries
%   xw - x-coordinates of water surface boundaries
%   yw - y-coordinates of water surface boundaries
%
%   * Patch coordinates (for use with patch plotting function)
%   n_patch - n-coordinate of edges for use with patch function
%   x_patch - x-coordinate of edges for use with patch function
%   y_patch - y-coordinate of edges for use with patch function
%   z_patch - z-coordinate of edges for use with patch function
%
%   * Indexing
%   col_to_mat - map column based data to matrix layout 
%   row_to_mat - map row based data to matrix layout
%   mat_to_cell - map matrix to cell vector layout
%   cell_to_mat - map cell vector layout to matrix layout
%   row_to_cell - map row based data to cell layout
%   col_to_cell - map column based data to cell layout
%
%   SigmaZetaMesh methods:
%   index - returns mesh cell indices for given positions
%   plot - plot the mesh optionally coloring with a given variable
%   plot3 - plot the mesh optionally coloring with a given variable in 3D
%
%   see also: Mesh, VelocitySolver


    properties (SetAccess=?SigmaZetaMeshGenerator)
        xs (1,1) XSection
        water_level (1,1) double 
        
        z_bottom_left (:,1) double {mustBeFinite, mustBeReal}
        z_top_left (:,1) double {mustBeFinite, mustBeReal}
        z_bottom_mid (:,1) double {mustBeFinite, mustBeReal}
        z_top_mid (:,1) double {mustBeFinite, mustBeReal}
        z_bottom_right (:,1) double {mustBeFinite, mustBeReal}
        z_top_right (:,1) double {mustBeFinite, mustBeReal}
        n_left (1,:) double {mustBeFinite, mustBeReal}
        n_middle (1,:) double {mustBeFinite, mustBeReal}
        n_right (1,:) double {mustBeFinite, mustBeReal}
        
        zb_left (1,:) double {mustBeFinite, mustBeReal}
        zb_middle (1,:) double {mustBeFinite, mustBeReal}
        zb_right (1,:) double {mustBeFinite, mustBeReal}
        zb_all (1,:) double {mustBeFinite, mustBeReal}
        nb_all (1,:) double {mustBeFinite, mustBeReal}
        
        nw (2,:) double {mustBeFinite, mustBeReal}
               
        col_to_mat (:,:) double {mustBeInteger, mustBeFinite mustBeReal}
        row_to_mat (:,:) double {mustBeInteger, mustBeFinite mustBeReal}
        mat_to_cell (:,1) logical
        cell_to_mat (:,1) double {mustBeInteger, mustBeFinite mustBeReal}
        row_to_cell (:,1) double {mustBeInteger, mustBeFinite mustBeReal}
        col_to_cell (:,1) double {mustBeInteger, mustBeFinite mustBeReal}
    end
    properties (Dependent, SetAccess=protected, GetAccess=public)
        x_left
        x_middle
        x_right
        y_left
        y_middle
        y_right
        xb_all
        yb_all
        xw
        yw
        n_patch
        x_patch
        y_patch
        z_patch
        sig_bottom_left
        sig_top_left
        sig_bottom_mid
        sig_top_mid
        sig_bottom_right
        sig_top_right
    end
    methods
        function val=get.x_left(obj)
            [val,~]=obj.xs.sn2xy(obj.n_left*0, obj.n_left);
        end
        function val=get.x_middle(obj)
            [val,~]=obj.xs.sn2xy(obj.n_middle*0, obj.n_middle);           
        end
        function val=get.x_right(obj)
            [val,~]=obj.xs.sn2xy(obj.n_right*0, obj.n_right);
        end
        function val=get.y_left(obj)
            [~,val]=obj.xs.sn2xy(obj.n_left*0, obj.n_left);
        end
        function val=get.y_middle(obj)
            [~,val]=obj.xs.sn2xy(obj.n_middle*0, obj.n_middle);           
        end
        function val=get.y_right(obj)
            [~,val]=obj.xs.sn2xy(obj.n_right*0, obj.n_right);
        end
        function val=get.xb_all(obj)
            nall=obj.nb_all;
            [val,~]=obj.xs.sn2xy(nall*0, nall);
        end
        function val=get.yb_all(obj)
            nall=obj.nb_all;
            [~,val]=obj.xs.sn2xy(nall*0, nall);
        end
        function val=get.xw(obj)
            nwl=obj.nw;
            [val,~]=obj.xs.sn2xy(nwl*0, nwl);
        end
        function val=get.yw(obj)
            nwl=obj.nw;
            [~,val]=obj.xs.sn2xy(nwl*0, nwl);
        end
        function val=get.x_patch(obj)
            npatch=obj.n_patch;
            [val,~]=obj.xs.sn2xy(npatch*0, npatch);
        end
        function val=get.y_patch(obj)
            npatch=obj.n_patch;
            [~,val]=obj.xs.sn2xy(npatch*0, npatch);
        end
        function val=get.sig_bottom_left(obj)
            val=obj.z_to_sigma(obj.z_bottom_left, reshape(obj.zb_left(obj.col_to_cell),[],1));
        end
        function val=get.sig_top_left(obj)
            val=obj.z_to_sigma(obj.z_top_left, reshape(obj.zb_left(obj.col_to_cell),[],1));
        end
        function val=get.sig_bottom_mid(obj)
            val=obj.z_to_sigma(obj.z_bottom_mid, reshape(obj.zb_middle(obj.col_to_cell),[],1));
        end
        function val=get.sig_top_mid(obj)
            val=obj.z_to_sigma(obj.z_top_mid, reshape(obj.zb_middle(obj.col_to_cell),[],1));
        end
        function val=get.sig_bottom_right(obj)
            val=obj.z_to_sigma(obj.z_bottom_right, reshape(obj.zb_right(obj.col_to_cell),[],1));
        end
        function val=get.sig_top_right(obj)
            val=obj.z_to_sigma(obj.z_top_right, reshape(obj.zb_right(obj.col_to_cell),[],1));
        end
        function s=z_to_sigma(obj,z, zb)
            s=(z-zb)./(obj.water_level-zb);
        end
        function val=get.n_patch(obj)
            val=[obj.n_left(obj.col_to_cell)
                 obj.n_middle(obj.col_to_cell)
                 obj.n_right(obj.col_to_cell)
                 obj.n_right(obj.col_to_cell)
                 obj.n_middle(obj.col_to_cell)
                 obj.n_left(obj.col_to_cell)
                 obj.n_left(obj.col_to_cell)];
        end
        function val=get.z_patch(obj)
            val=[obj.z_bottom_left,...
                 obj.z_bottom_mid,...
                 obj.z_bottom_right,...
                 obj.z_top_right,...
                 obj.z_top_mid,...
                 obj.z_top_left,...
                 obj.z_bottom_left]';
        end
        function idx=index(obj,n,sigma)
% Indices of mesh cells for given positions
%
%   idx = index(obj, n, sigma) returns the indices of the mesh cells that
%   hold the points given in (n,sigma) coordinates
%
%  see also: SigmaZetaMesh 
            idx=nan(size(sigma));
            nl=obj.n_left(obj.col_to_cell);
            nm=obj.n_middle(obj.col_to_cell);
            nr=obj.n_right(obj.col_to_cell);
            sbl=obj.sig_bottom_left;
            sbm=obj.sig_bottom_mid;
            sbr=obj.sig_bottom_right;
            stl=obj.sig_top_left;
            stm=obj.sig_top_mid;
            str=obj.sig_top_right;
            % TODO: vectorize code below!
            for cc=1:obj.ncells
                % left side
                fleft = n >= nl(cc) & n < nm(cc);
                fright = n >= nm(cc) & n < nr(cc);
                fsigleft = sigma > sbl(cc) + (sbm(cc) - sbl(cc))./ (nm(cc) - nl(cc)) .* n &...
                           sigma <= stl(cc) + (stm(cc) - stl(cc))./ (nm(cc) - nl(cc)) .* n &...
                           fleft;
                fsigright = sigma > sbm(cc) + (sbr(cc) - sbm(cc))./ (nr(cc) - nm(cc)) .* n &...
                            sigma <= stm(cc) + (str(cc) - str(cc))./ (nr(cc) - nm(cc)) .* n &...
                           fright;
                fincell = fsigleft | fsigright;
                idx(fincell) = cc;
            end
        end
        function plot(obj,var)
% Plot the mesh optionally colored with a variable
%
%   plot(obj) plots the mesh with the bed and water surface
%
%   plot(obj,var) plot the mesh and color the cells with the varibale var
%
%   see also: SigmaZetaMesh, plot3
            if ~isscalar(obj)
                for ce=1:numel(obj)
                    subplot(numel(obj),1,ce)
                    plot(obj(ce))
                end
                return
            end
            hold_stat=get(gca,'NextPlot');
            plot(obj.nb_all,obj.zb_all,'k','Linewidth',2)
            hold on
            plot(obj.nw,obj.nw*0+obj.water_level,'b','Linewidth',2)
            plot_var=nan(obj.ncells,1);
            if nargin > 1
                assert(numel(var)==obj.ncells, 'Variable to plot should have same number of elements as cells in the mesh');
                plot_var=var;
            end
            patch(obj.n_patch, obj.z_patch, plot_var(:));
            set(gca,'NextPlot',hold_stat);
        end
        function plot3(obj,var)
% Plot the mesh optionally colored with a variable in 3D
%
%   plot3(obj) plots the mesh with the bed and water surface
%
%   plot3(obj,var) plot the mesh and color the cells with the varibale var
%
%   see also: SigmaZetaMesh, plot
            if ~isscalar(obj)
                hold_stat=get(gca,'NextPlot');
                hold on
                for ce=1:numel(obj)
                    plot3(obj(ce))
                end
                set(gca,'NextPlot',hold_stat);
                return
            end
            hold_stat=get(gca,'NextPlot');
            plot3(obj.xb_all,obj.yb_all,obj.zb_all, 'k', 'Linewidth',2)
            hold on
            plot3(obj.xw,obj.yw,obj.water_level+obj.xw*0, 'b', 'Linewidth',2)
            plot_var=nan(obj.ncells,1);
            if nargin > 1
                assert(numel(var)==obj.ncells, 'Variable to plot should have same number of elements as cells in the mesh');
                plot_var=var;
            end
            patch(obj.x_patch,obj.y_patch, obj.z_patch,plot_var(:));
            set(gca,'NextPlot',hold_stat);
            set(gca,'DataAspectRatio',[5 5 1])
        end
        
    end
    methods(Access=protected)
        function val=get_ncells(obj)
            val=size(obj.z_bottom_left,1);
        end
    end
    
end