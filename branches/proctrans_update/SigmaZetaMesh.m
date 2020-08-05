classdef SigmaZetaMesh < Mesh
    properties
        % cross section
        xs (1,1) XSection
        water_level (1,1) double 
        
        % cell position
        z_bottom_left (:,1) double {mustBeFinite, mustBeReal}
        z_top_left (:,1) double {mustBeFinite, mustBeReal}
        z_bottom_mid (:,1) double {mustBeFinite, mustBeReal}
        z_top_mid (:,1) double {mustBeFinite, mustBeReal}
        z_bottom_right (:,1) double {mustBeFinite, mustBeReal}
        z_top_right (:,1) double {mustBeFinite, mustBeReal}
        n_left (1,:) double {mustBeFinite, mustBeReal}
        n_middle (1,:) double {mustBeFinite, mustBeReal}
        n_right (1,:) double {mustBeFinite, mustBeReal}
        
        % bed position
        zb_left (1,:) double {mustBeFinite, mustBeReal}
        zb_middle (1,:) double {mustBeFinite, mustBeReal}
        zb_right (1,:) double {mustBeFinite, mustBeReal}
               
        % indexing
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
        n_all
        z_all
        x_all
        y_all
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
        function val=get.x_all(obj)
            nall=obj.n_all;
            [val,~]=obj.xs.sn2xy(nall*0, nall);
        end
        function val=get.y_all(obj)
            nall=obj.n_all;
            [~,val]=obj.xs.sn2xy(nall*0, nall);
        end
        function val=get.x_patch(obj)
            npatch=obj.n_patch;
            [val,~]=obj.xs.sn2xy(npatch*0, npatch);
        end
        function val=get.y_patch(obj)
            npatch=obj.n_patch;
            [~,val]=obj.xs.sn2xy(npatch*0, npatch);
        end
        function val=get.n_all(obj)
            val=[reshape([obj.n_left;obj.n_middle],1,[]) obj.n_right(end)];
        end
        function val=get.z_all(obj)
            val=[reshape([obj.zb_left;obj.zb_middle],1,[]) obj.zb_right(end)];
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
            hold_stat=get(gca,'NextPlot');
            plot(obj.n_all,obj.z_all,'k','Linewidth',2)
            hold on
            plot_var=nan(obj.ncells,1);
            if nargin > 1
                assert(numel(var)==obj.ncells, 'Variable to plot should have same number of elements as cells in the mesh');
                plot_var=var;
            end
            patch(obj.n_patch, obj.z_patch, plot_var(:));
            set(gca,'NextPlot',hold_stat);
        end
        function plot3(obj,var)
            hold_stat=get(gca,'NextPlot');
            plot3(obj.x_all,obj.y_all,obj.z_all, 'k', 'Linewidth',2)
            hold on
            plot_var=nan(obj.ncells,1);
            if nargin > 1
                assert(numel(var)==obj.ncells, 'Variable to plot should have same number of elements as cells in the mesh');
                plot_var=var;
            end
            patch(obj.x_patch,obj.y_patch, obj.z_patch,plot_var(:));
            set(gca,'NextPlot',hold_stat);
            set(gca,'DataAspectRatio',[20 20 1])
        end
        
    end
    methods(Access=protected)
        function val=get_ncells(obj)
            val=size(obj.z_bottom_left,1);
        end
    end
    
end