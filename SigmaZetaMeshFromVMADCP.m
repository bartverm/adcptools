classdef SigmaZetaMeshFromVMADCP < SigmaZetaMeshGenerator
    properties
        vmadcp (1,1) VMADCP
        filter (1,1) EnsembleFilter
        bathymetry (1,1) Bathymetry = BathymetryScatteredPoints
        deltan (1,1) double = 5
        deltaz (1,1) double = 1
        xs (1,1) XSection
    end
    methods
        function obj=SigmaZetaMeshFromVMADCP(varargin)   
            construct_bathymetry=true;
            construct_xs=true;
            construct_filter=true;
            has_vmadcp=false;
            for count_arg=1:nargin         
                cur_arg=varargin{count_arg};
                if isa(cur_arg,'VMADCP')
                    has_vmadcp=true;
                    obj.vmadcp=cur_arg;
                elseif isa(cur_arg,'Bathymetry')
                    construct_bathymetry=false;
                    obj.bathymetry=cur_arg;
                elseif isa(cur_arg,'XSection')
                    construct_xs=false;
                    obj.xs=cur_arg;
                elseif isa(cur_arg,'EnsembleFilter')
                    construct_filter=false;
                    obj.filter=cur_arg;
                else
                    warning(['Unhandled input of type: ',class(cur_arg)])
                    % warning here
                end
            end
            if ~has_vmadcp
                error('Please pass a VMADCP object upon construction')
            end
            if construct_bathymetry
                obj.bathymetry=BathymetryScatteredPoints(obj.vmadcp);
            end
            if construct_xs
                obj.xs=XSection(obj.vmadcp);
            end
            if construct_filter
                obj.filter=EnsembleFilter(obj.vmadcp);
            end
        end
        function mesh=get_mesh(obj)
            mesh=SigmaZetaMesh;
            mesh.xs=obj.xs;
            % get velocity position and project on track
            vpos = obj.vmadcp.depth_cell_position;
            vpos(:,obj.filter.bad_ensembles,:,:) = nan;
            [~, vpos_n] = obj.xs.xy2sn(vpos(:, :, :, 1), vpos(:, :, :, 2) );

            % make n positions
            nmin = nanmin(vpos_n, [], 'all');  % compute minimum n coordinate of depth cells
            nmax = nanmax(vpos_n, [], 'all');  % compute maximum n coordinate of depth cells
            nvec = nmin : obj.deltan : nmax + obj.deltan;  % compute n position of cell edges
            nmid = (nvec(2 : end) + nvec(1 : end - 1) ) / 2; % compute n position of cell centers
            [xvec, yvec] = obj.xs.sn2xy(nvec*0, nvec); % compute xy positions of cell edges
            [xmid, ymid] = obj.xs.sn2xy(nmid*0, nmid); % compute xy positions of cell centers
            
            % interpolate bathymetry to n positions
            zvec = obj.bathymetry.get_bed_elev(xvec, yvec);                       % interpolate bathymetry to cell edges
            zmid = obj.bathymetry.get_bed_elev(xmid, ymid);                       % interpolate bathymetry to cell edges
                       
            % compute vertical limits of mesh
            maxz=nanmax(vpos(:,:,:,3),[],'all');
            minsigma=1-cosd(nanmax(obj.vmadcp.beam_angle));
            minz=zvec*(1-minsigma);
            minz_mid=zmid*(1-minsigma);
                       
            % create indices
            nz=max(0,ceil((maxz-minz_mid)/obj.deltaz));
            max_num=max(nz);
            mesh.col_to_mat=repmat(1:size(nz,2),max_num,1);
            mesh.row_to_mat=repmat((1:max_num)',1,size(nz,2));
            mesh.mat_to_cell=reshape(mesh.row_to_mat<=nz,[],1);
            mesh.row_to_cell=mesh.row_to_mat(mesh.mat_to_cell);
            mesh.col_to_cell=mesh.col_to_mat(mesh.mat_to_cell);
            mesh.cell_to_mat=reshape(1:numel(mesh.mat_to_cell),[],1);
            mesh.cell_to_mat(~mesh.mat_to_cell)=[];

            % store n-coordinates cells
            mesh.n_middle=reshape(nmid,1,[]);
            mesh.n_left=reshape(nvec(1:end-1),1,[]);
            mesh.n_right=reshape(nvec(2:end),1,[]);

            % store bed levels
            mesh.zb_middle=zmid;
            mesh.zb_left=zvec(1:end-1);
            mesh.zb_right=zvec(2:end);

            % create vertical position of cells
            minz_left=minz(1:end-1);
            minz_right=minz(2:end);
%             ncells=sum(nz);
            deltaz_mid=(maxz-minz_mid)./nz;
            deltaz_left=(maxz-minz_left)./nz;
            deltaz_right=(maxz-minz_right)./nz;
            mesh.z_bottom_left=reshape(maxz-mesh.row_to_mat.*deltaz_left,[],1);
            mesh.z_top_left=reshape(maxz-(mesh.row_to_mat-1).*deltaz_left,[],1);
            mesh.z_bottom_mid=reshape(maxz-mesh.row_to_mat.*deltaz_mid,[],1);
            mesh.z_top_mid=reshape(maxz-(mesh.row_to_mat-1).*deltaz_mid,[],1);
            mesh.z_bottom_right=reshape(maxz-mesh.row_to_mat.*deltaz_right,[],1);
            mesh.z_top_right=reshape(maxz-(mesh.row_to_mat-1).*deltaz_right,[],1);
            mesh.z_bottom_left=mesh.z_bottom_left(mesh.mat_to_cell);
            mesh.z_top_left=mesh.z_top_left(mesh.mat_to_cell);
            mesh.z_bottom_mid=mesh.z_bottom_mid(mesh.mat_to_cell);
            mesh.z_top_mid=mesh.z_top_mid(mesh.mat_to_cell);
            mesh.z_bottom_right=mesh.z_bottom_right(mesh.mat_to_cell);
            mesh.z_top_right=mesh.z_top_right(mesh.mat_to_cell);
            
        end
    end
end