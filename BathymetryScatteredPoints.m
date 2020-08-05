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
        function plot(obj)
           scatter3(obj.x,obj.y,obj.z,1,obj.z,'filled') 
           set(gca,'dataaspectratio',[20 20 1])
        end
    end
end