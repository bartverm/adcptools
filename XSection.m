classdef XSection < handle
    properties
        origin (2,1) double {mustBeFinite} = [0; 0]
        direction (2,1) double {mustBeFinite,XSection.mustBeUnitVector} = [1; 0];
        bed_position (2,:) double {mustBeFinite} = double.empty(2,0);
    end
    properties(Dependent)
        direction_orthogonal (2,1) double {mustBeFinite,XSection.mustBeUnitVector}
    end
    methods
        function val=get.direction_orthogonal(obj)
            val=[obj.direction(2); -obj.direction(1)];
        end
        function set.direction_orthogonal(obj,val)
            obj.direction=[-val(2) val(1)];
        end
        function [outpos,outvel]=xy2sn(obj,inpos,invel,dim)
            if (nargout == 2 && nargin < 3) || nargout > 2 
                error('Too many output arguments')
            end
            if nargin <4
                dim=1;
            end
            validateattributes(inpos,{'double'},{});
            assert(size(inpos,dim)==2,'Size of first dimension of inpos must be 2');
            if nargin > 2
                validateattributes(invel,{'double'},{});
                assert(size(invel,dim)==2,'Size of first dimension of invels must be 2');
            end
            outpos=cat(dim,XSection.dot(inpos-obj.origin,obj.direction_orthogonal,dim),...
                           XSection.dot(inpos-obj.origin,obj.direction,dim));
            if nargin > 2
                outvel=cat(dim,XSection.dot(invel,obj.direction_orthogonal,dim),...
                               XSection.dot(invel,obj.direction,dim));
            end
            
        end
        function [outpos,outvel]=sn2xy(obj,inpos,invel)
            if (nargout == 2 && nargin < 2) || nargout > 2 
                error('Too many output arguments')
            end
            validateattributes(inpos,{'double'},{'size',[2, nan]});
            validateattributes(invel,{'double'},{'size',[2, nan]});
            outpos=inpos(2,:).*obj.direction+obj.origin+inpos(1,:)*obj.direction_orthogonal;
            outvel=invel(2,:).*obj.direction+invel(1,:)*obj.direction_orthogonal;
        end
                
        function varargout=plot(obj,varargin)
            hbed=plot(obj.bed_position(1,:),obj.bed_position(2,:),varargin{:});
            if nargout>0
                varargout{1}=hbed;
            end
        end
    end
    methods (Static, Access=protected)
        function mustBeUnitVector(val)
            assert(sum(val.^2)-1<eps,'The value must be a unit vector')
        end
        function val=dot(a,b,dim)
            perm=1:max(ndims(a),ndims(b));
            perm(dim)=[];
            perm=[dim perm];
            val=ipermute(sum(permute(a,perm).*permute(b,perm)),perm);            
        end
    end
end