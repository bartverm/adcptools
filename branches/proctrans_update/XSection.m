classdef XSection < handle
    % Defines a cross-section
    %
    %   obj=XSection() Construct default cross-section
    %
    %   obj=XSection(...) pass different objects to initialize the object. This
    %   allows to pass a VMADCP object and/or an EnsembleFilter object which
    %   will be used to estimate the origin and direction of the cross-section
    %   based on the VMADCP track, filterd with the EnsembleFilter object. This
    %   can also be done later using the set_from_vmadcp method.
    %
    %   XSection properties:
    %   origin - origin coordinates of the cross-section
    %   direction - tangential direction of the cross-section
    %
    %   XSection properties (read only):
    %   direction_orthogonal - orthogonal direction of the cross-section
    %
    %   XSection methods:
    %   xy2sn - transform xy coordinates of points and vectors to sn coordinates
    %   sn2xy - transform sn coordinates of points and vectors to xy coordinates
    %   plot - plots the tangential and orthogonal vectors at the origin
    %   set_from_vmadcp - sets the origin and direction based on the vmadcp track
    %
    %   see also: VMADCP, EnsembleFilter
    properties
        % XSection/origin
        %
        % Projected coordinates of the origin of the cross-section. This
        % must be a two element column vector with finite values. Default
        % is [0; 0].
        %
        % see also: XSection, direction, direction_orthogonal
        origin (2,1) double {mustBeFinite} = [0; 0]
        
        % Xsection/direction
        %
        % Unit vector pointing in the tangential direction of the
        % cross-section. This must be a two element column vector with
        % finite values. Default is [1; 0]
        %
        % see also: XSection, origin, direction_orthogonal
        direction (2,1) double {mustBeFinite,XSection.mustBeUnitVector} = [1; 0];
    end
    properties(Dependent)
        %XSection/direction_orthogonal (read only)
        %
        % Unit vector orthogonal to the tangential direction.
        %
        % see also: XSection, origin, direction
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
                obj.set_from_vmadcp(vadcp, filter)
            end
        end
        function val=get.direction_orthogonal(obj)
            val=[obj.direction(2); -obj.direction(1)];
        end
        function set.direction_orthogonal(obj,val)
            obj.direction=[-val(2) val(1)];
        end
        function [s,n,us,un]=xy2sn(obj,x,y,u,v)
            % Transform points and vectors from xy to sn coordinates
            %
            %   [s, n] = xy2sn(obj,x,y) transforms the projected point coordinates
            %   given in (x,y) to (s,n) coordinates being the distance
            %   orthogonal and along the cross-section.
            %
            %   [s,n,us,un]=xy2sn(obj,x,y,u,v) allows to transform a vector
            %   located at (x,y) with component (u,v) both in projected
            %   coordinates to the corresponding (s,n) locations with (us,un)
            %   components across and along the cross-section respectively.
            %
            %   see also: XSection, sn2xy
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
            % Transform points and vectors from sn to xy coordinates
            %
            %   [x, y] = sn2xy(obj,s,n) transforms the point coordinates
            %   given in (s,n), i.e. accross and along coordinates to projected
            %   geographic coordinates in (x,y).
            %
            %   [x, y, u, v] = xy2sn(obj, s, n, us, un) allows to transform a vector
            %   located at (s,n) with component (us,un) both in accros and along
            %   coordinates to the corresponding (s,n) locations with (us,un)
            %   components across and along the cross-section respectively.
            %
            %   see also: XSection, sn2xy
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
            % Plot orthogonal and tangential vector at cross-section origin
            %
            %   plot(obj) plot the unit vectors orthogonal and tangential to the
            %   cross-section at the cross-section origin.
            %
            %   plot(obj,scale) manually set the scale of the vectors
            %
            %   h = plot(...) returns the handle to the quivers
            %
            %   see also: XSection
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
        function set_from_vmadcp(obj,V,filter)
% Set the direction and origin of the cross-section based on the adcp track
%
%   set_from_vmadcp(obj,V) sets the direction of the by calculating the 
%   largest eigenvector of the covariance matrix of the positions. The 
%   origin is determined by setting first the origin to the mean of the 
%   cross-section. Than the track is rotated to s,n coordinates and the 
%   origin is set to the mid-point of the n-range.
%
%   set_from_vmadcp(obj,V,filter) allows to specify an EnsembleFilter to 
%   exclude parts of the track
%
%   see also: XSection, EnsembleFilter, VMADCP
            [x,y]=V.xy;
            if nargin > 2 && ~isempty(filter)
                assert(isa(filter,'EnsembleFilter'),'filter should be of class EnsembleFilter')
                fbad=filter.all_cells_bad(V);
                x(fbad)=nan;
                y(fbad)=nan;
            end
            C=nancov([x;y]');
            [eigvec,eigval]=eig(C);
            [~, fprinc]=max(sum(eigval,2));
            obj.direction=eigvec(:,fprinc);
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
    end
end