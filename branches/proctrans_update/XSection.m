classdef XSection < handle
% Defines a cross-section for velocity data processing
    properties
        origin (2,1) double {mustBeFinite} = [0; 0]
        direction (2,1) double {mustBeFinite,XSection.mustBeUnitVector} = [1; 0];
    end
    properties(Dependent)
        direction_orthogonal (2,1) double {mustBeFinite,XSection.mustBeUnitVector}
    end
    methods
        function obj=XSection(varargin)
            filter=EnsembleFilter.empty(1,0);
            construct_from_vmadcp=false;
            for ca=1:nargin
                cur_arg=varargin{ca};
                if isa(cur_arg,'VMADCP')
                    construct_from_vmadcp=true;
                    vadcp=cur_arg;
                elseif isa(cur_arg,'EnsembleFilter')
                    filter=[filter, cur_arg]; %#ok<AGROW>
                end
            end
            if construct_from_vmadcp
                obj.set_direction_from_vmadcp(vadcp, filter)
                obj.set_origin_from_vmadcp(vadcp, filter)
            end
        end
        function val=get.direction_orthogonal(obj)
            val=[obj.direction(2); -obj.direction(1)];
        end
        function set.direction_orthogonal(obj,val)
            obj.direction=[-val(2) val(1)];
        end
        function [s,n,us,un]=xy2sn(obj,x,y,u,v)
            if nargout == 0
                return
            elseif nargin == 5
                has_vel=true;
            elseif nargin == 3
                has_vel=false;
            else 
                error('Function takes either 3 or 5 inputs')
            end
            
            if has_vel
                validateattributes(u,{'numeric'},{});
                validateattributes(v,{'numeric'},{});
                assert(isequal(size(u),size(v)),'size of u and v should match')
            end
            validateattributes(x,{'numeric'},{});
            validateattributes(y,{'numeric'},{});
            assert(isequal(size(x),size(y)),'size of x and y should match')
            s = (x - obj.origin(1)) * obj.direction_orthogonal(1) + (y - obj.origin(2)) * obj.direction_orthogonal(2);
            n = (x - obj.origin(1)) * obj.direction(1) + (y - obj.origin(2)) * obj.direction(2);
            if has_vel
                us = u * obj.direction_orthogonal(1) + v * obj.direction_orthogonal(2);
                un = u * obj.direction(1) + v * obj.direction(2);
            end
        end
        function [x, y, u, v]=sn2xy(obj, s, n, us, un)
            if nargout == 0
                return
            elseif nargin == 5
                has_vel=true;
            elseif nargin == 3
                has_vel=false;
            else 
                error('Function takes either 3 or 5 inputs')
            end
            
            if has_vel
                validateattributes(us,{'numeric'},{});
                validateattributes(un,{'numeric'},{});
                assert(isequal(size(us),size(un)),'size of u and v should match')
            end
            validateattributes(s,{'numeric'},{});
            validateattributes(n,{'numeric'},{});
            assert(isequal(size(s),size(n)),'size of x and y should match')
            x = obj.origin(1) + obj.direction(1) * n + obj.direction_orthogonal(1) * s;
            y = obj.origin(2) + obj.direction(2) * n + obj.direction_orthogonal(2) * s;
            if has_vel
                u = obj.direction(1) * un + obj.direction_orthogonal(1) * us;
                v = obj.direction(2) * un + obj.direction_orthogonal(2) * us;
            end
        end
                
        function varargout=plot(obj,scale)
            scalearg={'autoscale','off'};
            if nargin < 2
                scale=1;
                scalearg{2}='on';
            end
            hold_stat=get(gca,'nextplot');
            h(1)=quiver(obj.origin(1),obj.origin(2),obj.direction(1)*scale,obj.direction(2)*scale,'linewidth',2,'color','r',scalearg{:});
            hold on
            h(2)=quiver(obj.origin(1),obj.origin(2),obj.direction_orthogonal(1)*scale,obj.direction_orthogonal(2)*scale,'linewidth',2,'color','g',scalearg{:});
            asp=get(gca,'dataaspectratio');
            set(gca,'dataaspectratio',[max(asp(1:2))*[1 1] asp(3)])
            set(gca,'nextplot',hold_stat)
            if nargout>0
                varargout{1}=h;
            end
        end
        function set_direction_from_vmadcp(obj,V,filter)
            [x,y]=V.xy;
            fbad=filter.all_cells_bad(V);
            x(fbad)=nan;
            y(fbad)=nan;
            C=nancov([x;y]');
            [eigvec,eigval]=eig(C);
            [~, fprinc]=max(sum(eigval,2));
            obj.direction=eigvec(:,fprinc);       
        end
        function set_origin_from_vmadcp(obj,V,filter)
            [x,y]=V.xy;
            fbad=filter.all_cells_bad(V);
            x(fbad)=nan;
            y(fbad)=nan;
            obj.origin=[nanmean(x); nanmean(y)];
            [~ , n]=obj.xy2sn(x, y);
            n_orig=(nanmax(n)+nanmin(n))/2;
            [obj.origin(1), obj.origin(2)]=obj.sn2xy(0, n_orig);
        end
    end
    methods (Static, Access=protected)
        function mustBeUnitVector(val)
            assert(sum(val.^2)-1<=eps,'The value must be a unit vector')
        end
        function val=dot(a,b,dim)
            perm=1:max(ndims(a),ndims(b));
            perm(dim)=[];
            perm=[dim perm];
            val=ipermute(sum(permute(a,perm).*permute(b,perm)),perm);            
        end
    end
end