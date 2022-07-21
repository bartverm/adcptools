classdef SigmaZetaMeshFromVMADCP < SigmaZetaMeshGenerator & helpers.ArrayOnConstruct
% Generates a SigmaZetaMesh based on VMADCP data
%
%   obj=SigmaZetaFromaVMADCP() constructs a default object
%
%   obj=SigmaZetaFromVMADCP(...) allows to pass properties to construct the
%   object with. Based on the class the following properties are defined:
%   VMADCP - vmadcp property
%   EnsembleFilter - filter property
%   Bathymetry - bathymetry property
%   XSection - xs property
%   datetime - time property
%
%   SigmaZetaFromVMADCP properties:
%   vmadcp - VMADCP object to construct mesh with
%   filter - To exclude data in the mesh construction
%   bathymetry - Defines the bathymetry
%   deltan - horizontal (along the cross-section) mesh resolution
%   deltaz - vertical mesh resolution 
%   xs - Defines the cross-section
%   time - Defines the time for which mesh should be generated
%
%   SigmaZetaFromVMADCP methods:
%   get_mesh - construct the mesh
%
%   see also: SigmaZetaMesh, VMADCP
    properties
% SigmaZetaMeshFromVMADCP/vmadcp
%
%   scalar VMADCP object holding the adcp data
%
%   see also:SigmaZetaFromVMADCP, VMADCP
        vmadcp (1,1) VMADCP
        
% SigmaZetaMeshFromVMADCP/filter
%
%   scalar EnsembleFilter object to exclude ADCP data in the mesh
%   generation
%
%   see also:SigmaZetaFromVMADCP, EnsembleFilter
        filter (1,1) EnsembleFilter
        
% SigmaZetaMeshFromVMADCP/bathymetry
%
%   scalar Bathymetry object defining the bed elevation. Default value is
%   BathymetryScatteredPoints(vmadcp) which creates a bathymetry based on the
%   ADCP data.
%
%   see also:SigmaZetaFromVMADCP, Bathymetry
        bathymetry (1,1) Bathymetry = BathymetryScatteredPoints
        
% SigmaZetaMeshFromVMADCP/deltan
%
%   scalar, positive double defining the horizontal (along the
%   cross-section) resolution of the mesh. Default is 5.
%
%   see also:SigmaZetaFromVMADCP
        deltan (1,1) double {mustBeFinite, mustBePositive}= 5
        
% SigmaZetaMeshFromVMADCP/
%
%   scalar, positive double defining the vertical resolution of the mesh
%   in m. Default is 1.
%
%   see also:SigmaZetaFromVMADCP
        deltaz (1,1) double {mustBeFinite, mustBePositive}= 1
        
% SigmaZetaMeshFromVMADCP/xs
%
%   Scalar XSection object defining the cross-section. By default this is
%   XSection(vmadcp) which constructs a cross-section based on the vmadcp
%   track data.
%
%   see also:SigmaZetaFromVMADCP, XSection
        xs (1,1) XSection
        
% SigmaZetaMeshFromVMADCP/time
%
%   Scalar datetime object defining the time at which the mesh should be
%   generated. This determines the water level to be used. By default this
%   is the current time. This works fine for some ConstantWaterLevel, but
%   needs to be properly set for other types of WaterLevel.
%
%   see also: SigmaZetaFromVMADCP, WaterLevel, datetime
        time (1,1) datetime
    end
    properties (Dependent, SetAccess = protected)
% SigmaZetaMeshFromVMADCP/water_level
%
%   Read only property retrieving the water level at the given time.
%
%   see also: SigmaZetaFromVMADCP, datetime, WaterLevel
        water_level
    end
    methods
        function obj=SigmaZetaMeshFromVMADCP(varargin)
            obj = obj@helpers.ArrayOnConstruct(varargin{:});
            construct_bathymetry=true;
            construct_xs=true;
            construct_filter=true;
            construct_time=true;
            has_vmadcp=false;
            for count_arg=1:nargin         
                cur_arg=varargin{count_arg};
                if isa(cur_arg,'VMADCP')
                    has_vmadcp=true;
                    var_name='vmadcp';
                elseif isa(cur_arg,'Bathymetry')
                    construct_bathymetry=false;
                    var_name='bathymetry';
                elseif isa(cur_arg,'XSection')
                    construct_xs=false;
                    var_name='xs';
                elseif isa(cur_arg,'EnsembleFilter')
                    construct_filter=false;
                    var_name='filter';
                elseif isa(cur_arg,'datetime')
                    var_name='time';
                    construct_time=false;
                else
                    warning(['Unhandled input of type: ',class(cur_arg)])
                    continue
                    % warning here
                end
                obj.assign_var(var_name,cur_arg);
            end

            if ~has_vmadcp
                return
            end
            vm = [obj.vmadcp];

            cell_vm = num2cell(vm);
            scalar_vm = vm;
            if ~isscalar(cell_vm) && isequal(cell_vm{:})
                scalar_vm = vm(1);
            end
            
            if construct_bathymetry
                B=BathymetryScatteredPoints(scalar_vm);
                obj.assign_var('bathy',B);
            end
            if construct_xs
                obj.assign_var('xs',XSection(scalar_vm));
            end
            if construct_filter
                obj.assign_var('filter',EnsembleFilter(scalar_vm));
            end
            if construct_time
                t={vm.time};
                ef=[obj.filter];
                ef={ef.bad_ensembles};
                ef = cellfun(@not,ef,"UniformOutput",false);
                t = cellfun(@(a,b) mean(a(b),"all","omitnat"),t,ef,'UniformOutput',false);
                obj.assign_var('time',t);
            end


        end
        function val=get.water_level(obj)
            val=obj.vmadcp.water_level_object.get_water_level(obj.time);
        end
        function mesh=get_mesh(obj)
% Construct the SigmaZetaMesh
%
%   mesh = get_mesh(obj) construct a SigmaZetaMesh. The generated mesh
%   spans from the minimum to the maximum n coordinate in the horizontal
%   and vertically between a minimum sigma coordinate, where adcp side-lobe
%   interference is expected (determined based on the adcp beam angle) and
%   the highest depth cell in the data.
%
%   see also: SigmaZetaMeshFromVMADCP, SigmaZetaMesh, VMADCP

            % handle call from array
            if ~isscalar(obj)
                siz_obj = size(obj);
                siz_obj = num2cell(siz_obj);
                mesh(siz_obj{:})=SigmaZetaMesh;
                for co = 1:numel(obj)
                    mesh(co)=obj(co).get_mesh();
                end
                return
            end

            mesh=SigmaZetaMesh;
            mesh.xs=obj.xs;
            mesh.time=obj.time;

            % make n positions
            bpos=obj.vmadcp.bed_position();
            bpos(:,obj.filter.all_cells_bad(obj.vmadcp),:,:)=nan;
            [~,n]=obj.xs.xy2sn(bpos(:,:,:,1),bpos(:,:,:,2)); % compute n positions of bed detections
            [xn,yn]=obj.xs.sn2xy(n*0,n); % compute x,y coordinates of bed detections projected on xs line
            bpos=reshape(bpos,[prod(size(bpos,2,3)),size(bpos,4)]);
            as=alphaShape(bpos(all(isfinite(bpos),2),1:2),'HoleThreshold',obj.xs.scale.^2); % get alphashape for bed positions, removing small holes
            as.Alpha=as.criticalAlpha('one-region'); % set alpha radius to get one region
            fgood=isfinite(n);
            n=n(fgood);
            xn=xn(fgood);
            yn=yn(fgood);
            outside_as=~as.inShape(xn,yn);
            n(outside_as)=[]; % remove n-points outside alpha shape
            nmin = min(n, [], 'all','omitnan');  % compute minimum n coordinate of depth cells
            nmax = max(n, [], 'all','omitnan');  % compute maximum n coordinate of depth cells
            nvec = [nmin : obj.deltan : nmax nmax];  % compute n position of cell edges
            
            
            % interpolate bathymetry to n positions
            [xvec, yvec] = obj.xs.sn2xy(nvec*0, nvec); % compute xy positions of cell edges
            zvec = obj.bathymetry.get_bed_elev(xvec, yvec);                       % interpolate bathymetry to cell edges
            mesh.nb_all=nvec;
            mesh.zb_all=zvec;
            
            % compute vertical limits of mesh
            minsigma=1-cosd(max(obj.vmadcp.beam_angle,[],'omitnan'));
            maxz = sum(...
                interp1(obj.vmadcp.time', ...
                [ obj.vmadcp.vertical_position; ...
                  obj.vmadcp.distmidfirstcell; ...
                  obj.vmadcp.cellsize ]', ...
                obj.time) .* ...
                [1 -1 .5]...
                ); % i.e vertpos - distmidfirtcell + cellsize/2

            % Check for bathymetry crossing the waterlevel
            wl=obj.water_level;
            mesh.water_level=wl;
            fwl=SigmaZetaMeshFromVMADCP.get_intersections(zvec,wl); % find intersections of bed with water level
            if isempty(fwl)
                error('SigmaZetaMeshFromVMADCP:get_mesh:WaterLevelBelowBed','Cannot create mesh since the water level is lower than the bed')
            end
            mesh.nw=reshape(nvec(fwl),2,[]);
            fi= fwl~=1  & fwl~=(numel(nvec)) | ( fwl==1 & reshape(zvec(fwl), 2,[]) > wl ); % find cells for which n-position of wl needs fixing. This is never the last point (so fwl+1 is fine below)
            fwl=fwl(fi);
            mesh.nw(fi)=nvec(fwl) + (nvec(fwl+1) - nvec(fwl)) ./...
                                        (zvec(fwl+1) - zvec(fwl)) .*...
                                        (wl - zvec(fwl));
            
            % Check for lowest mesh elevation crossing max mesh level
            z_bot_mesh=zvec*(1-minsigma)+minsigma*wl;
            
            fbnds=SigmaZetaMeshFromVMADCP.get_intersections(z_bot_mesh,maxz); % find intersections between minimum measurement level and maximum mesh level
            if isempty(fbnds)
                error('SigmaZetaMeshFromVMADCP:get_mesh:MaxZBelowBed','Cannot create mesh since the maximum mesh level is lower than the minimum measurement level')
            end
            
            % fix degenerate cells crossing maximum mesh level
            n_left=reshape(nvec(1:end-1),1,[]); % create n-coordinate of left side
            n_right=reshape(nvec(2:end),1,[]); % create n-coordinate of right side
            lf=fbnds(1,fbnds(1,:) ~= 1 | (fbnds(1,:) == 1 &  z_bot_mesh(fbnds(1,:))>maxz)); %find cells that need fixing on the left side
            n_left(lf)= n_left(lf) + (n_left(lf+1) - n_left(lf)) ./...
                                     (z_bot_mesh(lf+1) - z_bot_mesh(lf)) .*...
                                     (maxz-z_bot_mesh(lf)); % fix n_coordinate to match intersection point
            rf=fbnds(2,fbnds(2,:) <= numel(n_right) | (fbnds(1,:) == numel(n_right) & z_bot_mesh(fbnds(2,:))>maxz)); %find cells that need fixing on the right side
            n_right(rf)= n_right(rf-1) + (n_right(rf) - n_right(rf-1)) ./...
                                        (z_bot_mesh(rf) - z_bot_mesh(rf-1)) .*...
                                        (maxz-z_bot_mesh(rf)); % fix n_coordinate to match intersection point
            
            % compute bed levels at cell boundaries
            zb_left=zvec(1:end-1);
            zb_right=zvec(2:end);
            zb_left(lf) = (maxz - wl*minsigma)/(1-minsigma);
            zb_right(rf) = (maxz - wl*minsigma)/(1-minsigma);
                                    
            % Remove cells that are entirely above maximum mesh level
            frem=z_bot_mesh>maxz; % find boundaries where bottom is above maximum level
            frem=frem(1:end-1) & frem(2:end); % find cells where both left and right boundaries are above maxz
            n_left(frem)=[]; 
            n_right(frem)=[];
            zb_left(frem)=[];
            zb_right(frem)=[];
            
            % find bed elevation at the center of verticals
            nmid=(n_left+n_right)/2;
            [xmid, ymid] = obj.xs.sn2xy(nmid*0, nmid); % compute xy positions of cell centers            
            zmid = obj.bathymetry.get_bed_elev(xmid, ymid);                       % interpolate bathymetry to cell edges
            minz_mid=zmid*(1-minsigma)+minsigma*wl;
            [nb_all, idx_unq]=unique([mesh.nb_all, reshape(mesh.nw,1,[]), n_left, n_right, nmid]); % include on the bed line all points used on the bed. only unique, sorted points are inlcuded
            zb_all=[mesh.zb_all, reshape(mesh.nw,1,[])*0+wl, zb_left, zb_right, zmid];
            zb_all=zb_all(idx_unq);
            mesh.nb_all=nb_all;
            mesh.zb_all=zb_all;
                                 
            % create indices
            nz=max(1,ceil((maxz-minz_mid)/obj.deltaz)); % compute number of cells in each vertical
            max_num=max(nz); % get maximum number of cells in a vertical
            mesh.col_to_mat=repmat(1:size(nz,2),max_num,1);
            mesh.row_to_mat=repmat((1:max_num)',1,size(nz,2));
            mesh.mat_to_cell=reshape(mesh.row_to_mat<=nz,[],1);
            mesh.row_to_cell=mesh.row_to_mat(mesh.mat_to_cell);
            mesh.col_to_cell=mesh.col_to_mat(mesh.mat_to_cell);
            mesh.cell_to_mat=reshape(1:numel(mesh.mat_to_cell),[],1);
            mesh.cell_to_mat(~mesh.mat_to_cell)=[];

            % store n-coordinates cells
            mesh.n_middle=reshape(nmid,1,[]);
            mesh.n_left=reshape(n_left,1,[]);
            mesh.n_right=reshape(n_right,1,[]);

            % store bed levels
            mesh.zb_middle=zmid;
            mesh.zb_left=zb_left;
            mesh.zb_right=zb_right;

            % create vertical position of cells
            minz_left=zb_left*(1-minsigma)+minsigma*wl;
            minz_right=zb_right*(1-minsigma)+minsigma*wl;
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
    methods(Access=protected, Static)
        function f=get_intersections(vec,lev)
            if vec(1) < lev
                fstart=1;
            else
                
                fstart=[];
            end
            if vec(end)<lev
                fend=numel(vec);
            else
                fend=[];
            end
            fstart=[fstart find(diff(vec > lev) == -1)];
            if isempty(fstart)
                f=double.empty(2,0);
                return
            end
            fend=[find(diff(vec > lev) == 1) fend];
            f=[fstart; fend];
        end
    end
end