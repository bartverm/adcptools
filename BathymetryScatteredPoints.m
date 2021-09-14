classdef BathymetryScatteredPoints < Bathymetry
    properties
        x (1,:) double {mustBeFinite} = double.empty(1,0);
        y (1,:) double {mustBeFinite} = double.empty(1,0);
        z (1,:) double {mustBeFinite} = double.empty(1,0);
        interpolator (1,1) Interpolator = LoessInterpolator();
    end
    methods
        function obj=BathymetryScatteredPoints(varargin)
            obj=obj@Bathymetry(varargin{:});
            filter=EnsembleFilter.empty(1,0);
            construct_from_vmadcp=false;
            for ca=1:nargin
                cur_arg=varargin{ca};
                if isa(cur_arg,'VMADCP')
                    construct_from_vmadcp=true;
                    vadcp=cur_arg;
                elseif isa(cur_arg,'EnsembleFilter')
                    filter=[filter, cur_arg]; %#ok<AGROW>
                elseif isa(cur_arg,'Interpolator')
                    obj.interpolator=cur_arg;
                else
                    warning('Bathymetry:unhadled_input',['Unhandled input of type: ', class(cur_arg)])
                end
            end
            if construct_from_vmadcp
                obj.pos_from_vmadcp(vadcp,filter)
            end
        end
        function z=get_bed_elev(obj,x,y)
            validateattributes(x,{'double'},{})
            validateattributes(y,{'double'},{})
            assert(isequal(size(x),size(y)),'Size of x and y should match')
            assert(isequal(size(obj.x),size(obj.y), size(obj.z)),'Size of x, y and z given in object should match')
            sizin=size(x);
            obj.interpolator.known = [obj.x; obj.y; obj.z];
            z=reshape(obj.interpolator.interpolate([reshape(x,1,[]); reshape(y,1,[])]),sizin);
        end
        function pos_from_vmadcp(obj,vmadcp,filter)
            validateattributes(vmadcp,{'VMADCP'},{'scalar'},'pos_from_vmadcp','vmadcp',2)
            tpos=vmadcp.bed_position;
            tpos(:,filter.bad_ensembles,:)=nan;
            xpos=tpos(:,:,:,1);
            ypos=tpos(:,:,:,2);
            zpos=tpos(:,:,:,3);
            isfin=all(isfinite(tpos),4);
            obj.x=xpos(isfin);
            obj.y=ypos(isfin);
            obj.z=zpos(isfin);
        end
        function plot_residuals(obj)
            z_interp=obj.get_bed_elev(obj.x,obj.y);
            res=z_interp-obj.z;
            scatter(obj.x,obj.y,5,res,'filled')
            set(gca,'clim',nanmean(res) + nanstd(res)*2*[-1 1])
            colorbar
            axis equal
        end
        function varargout=plot(obj)
           if ~isscalar(obj)
               plot@Bathymetry(obj);
               return
           end
           hold_stat=get(gca,'NextPlot');
           hp=plot3(obj.x,obj.y,obj.z,'k.');
           hold on
           set(gca,'dataaspectratio',[5 5 1])
           as=alphaShape(obj.x',obj.y',1);
           as.Alpha=as.criticalAlpha('one-region')*6;
           tri = alphaTriangulation(as);
           z_interp=obj.get_bed_elev(obj.x,obj.y);
           ht=trimesh(tri,obj.x, obj.y, z_interp ,'FaceColor','interp','EdgeColor','none');
           set(gca,'NextPlot',hold_stat)
           if nargout > 0
               varargout{1}=hp;
           end
           if nargout > 1
               varargout{2}=ht;
           end
        end
    end
end