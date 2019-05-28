classdef MeshCell < handle
    properties
        coordinates(2,:) double {mustBeFinite} = double.empty(2,0)
    end
    properties(Dependent,SetAccess=protected, GetAccess=public)
        n_coordinates(1,1) double {mustBeInteger}
    end
    methods
        function obj=MeshCell(coordinates)
            if nargin > 0
                obj.coordinates=coordinates;
            end
        end
        function val=get.n_coordinates(obj)
            val=size(obj.coordinates,2);
        end
        function set.coordinates(obj,val)
            if ~isequal(val(:,1),val(:,end))
                val=[val val(:,1)];
            end
            obj.coordinates=val;
        end
        function tf=isInCell(obj,qcoord,dim)
            if isempty(obj)
                tf=[];
                return
            end
            if iscalar(obj)
                if nargin < 3
                    dim = 1;
                end
                assert(size(qcoord,dim)==2,['qcoord must have size 2 along dimension ', dim])
                perm=1:ndims(qcoord);
                perm(dim)=[];
                perm=[dim perm];
                qcoord=permute(qcoord,perm);
                siz_perm=size(qcoord);
                qcoord=reshape(qcoord,2,[]);
                tf=ipermute(reshape(inpolygon(qcoord(1,:),qcoord(2,:),obj.coordinates(1,:),obj.coordinates(2,:)),...
                    siz_perm), perm);
            else
                tf=nan(size(qcoord));
                for co=1:numel(obj)
                    tf(obj.isInCell(qcoord,dim))=co;
                end
            end
        end
        function hout=plot(obj,varargin)
            if isempty(obj) 
                hp=[];
            elseif isscalar(obj)
                patch(obj.coordinates(1,:),obj.coordinates(2,:),zeros(1,obj.n_coordinates),varargin{:})
            else
                hp=nan(size(obj));
                for co=1:numel(obj)
                    hp(co)=obj(co).plot(varargin{:});
                end
            end
            if nargout>0
                hout=hp;
            end
        end
        
    end
end