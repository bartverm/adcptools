classdef RotationMatrix < handle
    % RotationMatrix represents a 3D rotation matrix
    methods
        function obj=RotationMatrix(varargin)
            % RotationMatrix create a RotationMatrix object
            %    obj=RotationMatrix() creates a Rotation matrix object with an identity
            %       rotation matrix. Use obj(N,M,...) = RotationMatrix() to create an
            %       M x N x ... array of RotationMatrix objects.
            %
            %   obj=RotationMatrix(R) creates a RotationMatrix based on the array R. R
            %       is an MxNx...x3x3 array which creates an array of MxNx... objects.
            %       If R is a 3x3 array it will create one object. In case the size of
            %       R is Nx3x3 a 1xN array is created.
            %       If R is a RotationMatrix object it will create copies of the
            %       objects in R. The created objects do not share handles with the
            %       R. This is the same as obj=RotationMatrix(R.get())
            %
            %   obj=RotationMatrix(R, SIZE) allows to specify the size of the output.
            %       This is only usefull if R is a 3x3 matrix which is used for the
            %       construction of an array of RotationMatrix objects with size SIZE.
            %       In case R is MxNx...x3x3 construction will only succeed if SIZE is
            %       equal to [M N ...]. Rotation matrices in R should be orthogonal.
            %
            %   obj=RotationMatrix(ANGLE1, ANGLE2, ANGLE3) Construct an object given
            %       three subsequent rotation angles in radians. The size of the output
            %       array is equal to the size of the inputs. It is assumed that the
            %       angles describe a body rotation around the x, y and z axes,
            %       respectively.
            %
            %   obj=RotationMatrix(ANGLE1, ANGLE2, ANGLE3, ORDER) Allows to specify the
            %       order of rotation with four characters. The first three specify the
            %       axes order (e.g. 'xyz' or 'xyx') and are on of 'x', 'y' or 'z'. The
            %       last character is a 'b' or a 'w' describing a body (or vector)
            %       rotation or a world (or axis) rotation. Examples of valid inputs
            %       are 'xyzw', 'xyzb', 'xzxb'. If order is empty it will default to
            %       'xyzb'.
            %
            %   obj=RotationMatrix(ANGLE1, ANGLE2, ANGLE3, ORDER, UNIT) allows to
            %       specify in what unit the angles are given. UNIT can be either
            %       'radians', the default, or 'degrees'.
            if nargin == 1
                siz=size(varargin{1});
                if ~isa(varargin{1},'RotationMatrix')
                    siz=siz(1:end-2);
                end
                siz=num2cell(siz);
                obj(siz{:})=RotationMatrix;
            elseif nargin == 2
                siz=num2cell(varargin{2});
                obj(siz{:})=RotationMatrix;
            elseif nargin >= 3
                obj(size(varargin{1}))=RotationMatrix;
            end
            obj.set(varargin{:});
        end % RotationMatrix()
        
        function val=get(obj)
            % get returns the object (array) matrix(ces)
            %
            %   M=obj.get() returns the rotation matrix for the object or the matrices
            %       if obj is an array of objects. If obj is a single object a 3x3
            %       matrix is returned otherwise an array of size [size(obj) 3 3] is
            %       returned
            if isscalar(obj)
                val=obj.m_matrix;
            else
                val=reshape(permute(cat(3,obj.m_matrix),[3,1,2]),[size(obj), 3, 3]);
            end
        end % get()
        
        function set(obj,varargin)
            % set Sets the rotation matrix of a RotationMatrix object (array)
            %
            %   This function takes the exact same arguments as the object constructor
            if nargin == 1, return, end
            if nargin == 2 % Rotation matrices are given
                if isa(varargin{1},'RotationMatrix') % copy construction
                    rot=varargin{1}.get();
                else
                    rot=varargin{1};
                end
            elseif nargin == 3 % Rotation matrix and size is given
                rot=varargin{1};
            elseif nargin > 3
                [angle1, angle2, angle3]=varargin{1:3};
                if nargin>=5 && ~isempty(varargin{4})
                    order=varargin{4};
                else
                    order='xyzb';
                end
                if nargin==6
                    if ischar(varargin{5}) && strcmpi(varargin{5},'degrees')
                        angle1=deg2rad(angle1);
                        angle2=deg2rad(angle2);
                        angle3=deg2rad(angle3);
                    elseif ~(ischar(varargin{5}) && strcmpi(varargin{5},'radians'))
                        throw(MException('RotationMatrix:WrongInput5','Did not understand fifth input'))
                    end
                end
                rot=RotationMatrix.three_angle_rotation_(angle1,angle2,angle3,order);
            end % if nargin
            
            numel_obj=numel(obj);
            dims_rot=ndims(rot);
            if dims_rot == 2
                rot=shiftdim(rot,-1);
                rot=repmat(rot,[numel_obj,1]);
            else
                if dims_rot == 3
                    rot=permute(rot,[4,1,2,3]);
                end
                rot=reshape(rot,[numel_obj,3,3]);
            end
            assert(all(reshape(RotationMatrix.is_orthogonal_(rot),[],1)),'RotationMatrix:NonOrthogonal','Input matrix should be orthogonal')
            rot=reshape(rot,[],9);
            for cmat=1:size(rot,1)
                obj(cmat).m_matrix(:)=rot(cmat,:);
            end
        end % set()
        
        function [angle1, angle2, angle3]=decompose_rotation(obj,order,unit)
            % decompose_rotations decomposes the rotation in three subsequent angles
            %
            %   [angle1, angle2, angle3]=decompose_rotation(OBJ, ORDER) Decomposes the
            %       rotation matrices given in OBJ, into three subsequent rotations
            %       around the given axes. The three outputs ANGLE1, ANGLE2 and ANGLE3
            %       have the same size as OBJ. The angles are given in radians.
            %       ORDER allows to specify the order of the rotation angles. The first
            %       three characters specify the axes of rotation and the last
            %       character defines whether it is a world or a body rotation.
            %       The order of the axes matters. There are two types of
            %       decompsitions: 2 and 3 axes. This leads to the following 24
            %       possible values for ORDER:
            %         Three axes:
            %            World rotations:
            %               xyzw yzxw zxyw (positive axes ordering)
            %               xzyw yxzw zyxw (negative axes ordering)
            %            Body rotations:
            %               xyzb yzxb zxyb (positive axes ordering)
            %               xzyb yxzb zyxb (negative axes ordering)
            %         Two axes:
            %            World rotations:
            %               xyxw yzyw zxzw (positive axes ordering)
            %               xyxw yxyw zyzw (negative axes ordering)
            %            Body rotations:
            %               xyxb yzyb zxzb (positive axes ordering)
            %               xzxb yxyb zyzb (negative axes ordering)
            %
            %  [angle1, angle2, angle3]=decompose_rotation(OBJ, ORDER, UNIT) Allows to
            %      specify the angle unit.
            %
            %   Note that the first angle, represents the first rotation
            %   applied. If, for instance, we have an xyz rotation we would
            %   have the following order of matrix multiplication RZ*RY*RX,
            %   such that the RX transformation is applied first.
            assert(ischar(order) && numel(order)==4 && isrow(order),'RotationMatrix:WrongOrder','Order must be a 4 character row vector')
            assert(any(strcmp(order(1),{'x','y','z'})),'RotationMatrix:WrongOrder','First character must be ''x'', ''y'' or ''z''')
            assert(any(strcmp(order(2),{'x','y','z'})),'RotationMatrix:WrongOrder','Second character must be ''x'', ''y'' or ''z''')
            assert(any(strcmp(order(3),{'x','y','z'})),'RotationMatrix:WrongOrder','Fourth character must be ''x'', ''y'' or ''z''')
            assert(any(strcmp(order(4),{'w','b'})),'RotationMatrix:WrongOrder','Fourth character must be ''w'' or ''b''')
            
            pos=true;
            switch order(3)
                case 'x', ax=1;
                    switch order(2)
                        case 'x', throw(MException('RotationMatrix:WrongOrder','Second axis cannot be the same as the first'));
                        case 'z', pos=false;
                    end
                case 'y', ax=2;                                        % yz
                    switch order(2)
                        case 'y', throw(MException('RotationMatrix:WrongOrder','Second axis cannot be the same as the first'));
                        case 'x', pos=false;
                    end
                case 'z', ax=3;                                        % zx
                    switch order(2)
                        case 'z', throw(MException('RotationMatrix:WrongOrder','Second axis cannot be the same as the first'));
                        case 'y', pos=false;
                    end
            end
            if numel(intersect(order,'xyz'))==3
                [angle3, angle2, angle1]=helpers.RotationMatrix.decompose_three_axes_(ax,pos,obj.get());
            elseif numel(intersect(order,'xyz'))==2
                assert(strcmp(order(1),order(3)),'RotationMatrix:WrongOrder','In a two-axes decomposition the first and last axis must be the same')
                [angle3, angle2, angle1]=helpers.RotationMatrix.decompose_two_axes_(ax,pos,obj.get());
            else
                throw(MException('RotationMatrix:WrongOrder','Axes ordering is wrong. See help for valid axes orders'));
            end
            sgn=-RotationMatrix.wb2sign_(order(4));
            angle1=sgn*angle1;
            angle2=sgn*angle2;
            angle3=sgn*angle3;
            if nargin==3
                if strcmp(unit,'degrees')
                    angle1=rad2deg(angle1);
                    angle2=rad2deg(angle2);
                    angle3=rad2deg(angle3);
                elseif ~strcmp(unit,'radians')
                    throw(MException('RotationMatrix:WrongUnit','Did not understand the angle unit'))
                end
            end
        end
    end
    
    properties(Access=private)
        m_matrix=eye(3);
    end
    
    methods(Static)
        function [a1,a2,a3,pos]=translate_order_(order)
            assert(ischar(order),'RotationMatrix:OrderNotChar','The rotation order should be given as characters')
            assert(isvector(order),'RotationMatrix:OrderNotVector','The rotation order characters should be given as vector')
            assert(numel(order)==4,'RotationMatrix:OrderNot4Elements','The rotation order should be exactly for characters long')
            a1=helpers.RotationMatrix.axis2num_(order(1));
            a2=helpers.RotationMatrix.axis2num_(order(2));
            a3=helpers.RotationMatrix.axis2num_(order(3));
            pos=helpers.RotationMatrix.wb2sign_(order(4));
        end
        function sign=wb2sign_(wb)
            switch wb
                case 'w', sign=1;
                case 'b', sign=-1;
                otherwise, throw(MException('RotationMatrix:WbNotWb','Rotation reference should be ''w'' (world rotation) or ''b'' (body rotation)'))
            end
        end
        function num=axis2num_(axis)
            switch axis
                case 'x', num=1;
                case 'y', num=2;
                case 'z', num=3;
                otherwise, throw(MException('RotationMatrix:AxisNotXyz','Axis should be ''x'', ''y'' or ''z'''))
            end
        end
        function rot = generate_matrix_(angle,axis)
            c=cos(angle);
            s=sin(angle);
            if isscalar(angle) || isvector(angle)
                angle=angle(:);
                nd=3;
            else
                nd=ndims(angle)+2;
            end
            siz_in=size(angle);
            switch axis
                case 'x'
                    rot=cat(nd-1, cat(nd,  ones(siz_in), zeros(siz_in), zeros(siz_in) ),...
                        cat(nd, zeros(siz_in),             c,            -s ),...
                        cat(nd, zeros(siz_in),             s,             c ));
                case 'y'
                    rot=cat(nd-1, cat(nd,             c, zeros(siz_in),            s  ),...
                        cat(nd, zeros(siz_in),  ones(siz_in), zeros(siz_in) ),...
                        cat(nd,            -s, zeros(siz_in),             c ));
                case 'z'
                    rot=cat(nd-1, cat(nd,             c,            -s, zeros(siz_in) ),...
                        cat(nd,             s,             c, zeros(siz_in) ),...
                        cat(nd, zeros(siz_in), zeros(siz_in),  ones(siz_in) ));
            end
        end
        function rot = three_angle_rotation_(angle1,angle2,angle3,order)
            assert(isnumeric(angle1) && isnumeric(angle2) && isnumeric(angle3),'RotationMatrix:AnglesNotNumeric','Angles should be numeric')
            assert(isequal(size(angle1),size(angle2),size(angle3)),'RotationMatrix:AnglesDifferentSize','Angle arrays should have the same size')
            assert(all(isfinite(angle1(:))) && all(isfinite(angle2(:))) && all(isfinite(angle3(:))),'RotationMatrix:AnglesNotFinite','Angles should be finite')
            assert(ischar(order) && numel(order)==4,'RotationMatrix:OrderNot4Char', 'Order must be given as a four character variable')
            siz_out=size(angle1);
            angle1=reshape(angle1,[],1);
            angle2=reshape(angle2,[],1);
            angle3=reshape(angle3,[],1);
            R1=helpers.RotationMatrix.generate_matrix_(helpers.RotationMatrix.wb2sign_(order(4))*angle1,order(1));
            R2=helpers.RotationMatrix.generate_matrix_(helpers.RotationMatrix.wb2sign_(order(4))*angle2,order(2));
            R3=helpers.RotationMatrix.generate_matrix_(helpers.RotationMatrix.wb2sign_(order(4))*angle3,order(3));
            rot=helpers.RotationMatrix.mult_(R3,helpers.RotationMatrix.mult_(R2,R1));
            rot=reshape(rot,[siz_out, 3,3]);
        end
        function b=transp_(a)
            b=permute(a,[1 3 2]);
        end
        function c=mult_(a,b)
            c=zeros([size(a,1),size(a,2),size(b,3)]);
            for cr=1:size(a,2)
                for cc=1:size(b,3)
                    for ce=1:size(a,3)
                        c(:,cr,cc)=c(:,cr,cc)+a(:,cr,ce).*b(:,ce,cc);
                    end
                end
            end
        end
        function yn=is_orthogonal_(rot)
            prod=helpers.RotationMatrix.mult_(rot,helpers.RotationMatrix.transp_(rot));
            yn=all(reshape((prod-repmat(shiftdim(eye(3),-1),[size(prod,1) 1])).^2<eps,[],1));
        end
        function [angle1, angle2, angle3]=decompose_three_axes_(axis_1, pos, M)
            % Herter and Lott 1993, Algorithm for decomposing 3-D orthogonal matrices
            % hinto primitive rotations.
            if ismatrix(M)
                shiftdim(M,-1);
                siz_out=[1,1];
            else
                siz_out=size(M);
                siz_out=siz_out(1:end-2);
            end
            
            if pos
                axis_2=mod(axis_1,3)+1;
                axis_3=mod(axis_2,3)+1;
            else
                axis_2=mod(axis_1-2,3)+1;
                axis_3=mod(axis_2-2,3)+1;
            end
            
            M=reshape(M,[],3,3);
            angle2=asin(M(:,axis_1,axis_3));
            if pos
                angle2=-angle2;
            end
            
            [angle3, sin_phi, cos_phi]=deal(nan(size(angle2)));
            is_sing= abs(M(:,axis_1,axis_3))<(1.0-eps);
            sin_phi(is_sing)=M(is_sing,axis_1,axis_2);
            if ~pos
                sin_phi(is_sing)=-sin_phi(is_sing);
            end
            angle3(is_sing)=atan2(reshape(sin_phi(is_sing),sum(is_sing),1),reshape(M(is_sing,axis_1,axis_1),sum(is_sing),1));
            sin_phi(is_sing)=-M(is_sing,axis_2,axis_3);
            cos_phi(is_sing)=M(is_sing,axis_3,axis_3);
            angle3(~is_sing)=0;
            sin_phi(~is_sing)=M(~is_sing,axis_3,axis_2);
            cos_phi(~is_sing)=M(~is_sing,axis_2,axis_2);
            if pos
                sin_phi=-sin_phi;
            end
            angle1=atan2(sin_phi,cos_phi);
            
            angle1=reshape(angle1,siz_out);
            angle2=reshape(angle2,siz_out);
            angle3=reshape(angle3,siz_out);
        end % decompos_3_axes
        
        function [angle1, angle2, angle3]=decompose_two_axes_(axis_1, pos, M)
            % Herter and Lott 1993, Algorithm for decomposing 3-D orthogonal matrices
            % hinto primitive rotations.
            if ismatrix(M)
                M=shiftdim(M,-1);
                siz_out=[1,1];
            else
                siz_out=size(M);
                siz_out=siz_out(1:end-2);
            end
            
            if pos
                axis_2=mod(axis_1,3)+1;
                axis_3=mod(axis_2,3)+1;
            else
                axis_2=mod(axis_1-2,3)+1;
                axis_3=mod(axis_2-2,3)+1;
            end
            
            M=reshape(M,[],3,3);
            angle2=acos(M(:,axis_1,axis_1));
            
            [angle3, sin_phi, cos_phi]=deal(nan(size(angle2)));
            is_sing= abs(M(:,axis_1,axis_1))<(1.0-eps);
            cos_phi(is_sing)=M(is_sing,axis_1,axis_3);
            if pos
                cos_phi(is_sing)=-cos_phi(is_sing);
            end
            angle3(is_sing)=atan2(reshape(M(is_sing,axis_1,axis_2),[sum(is_sing),1]),reshape(cos_phi(is_sing),[sum(is_sing),1]));
            sin_phi(is_sing)=M(is_sing,axis_2,axis_1);
            cos_phi(is_sing)=M(is_sing,axis_3,axis_1);
            if ~pos
                cos_phi(is_sing)=-cos_phi(is_sing);
            end
            angle3(~is_sing)=0;
            cos_phi(~is_sing)=M(~is_sing,axis_2,axis_2);
            sin_phi(~is_sing)=M(~is_sing,axis_3,axis_2);
            if pos
                sin_phi(~is_sing)=-sin_phi(~is_sing);
            end
            angle1=atan2(sin_phi,cos_phi);
            
            angle1=reshape(angle1,siz_out);
            angle2=reshape(angle2,siz_out);
            angle3=reshape(angle3,siz_out);
        end % decompos_two_axes
    end
end