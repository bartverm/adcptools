function out=matmult(a,b,dim1,dim2)
% matmult performs matrix multiplications on multidimensional arrays
%
%   OUT=matmult(A,B) performs matrix multiplication over A and B assuming
%   the matrices are defined in the last two dimensions. Inner dimensions
%   must be equal. All dimensions other than those in which the matrices
%   are defined have to be equally sized or singleton.
%   In this case: DIM2=max(ndims(A),ndims(B)) and DIM1=DIM2-1;
%   
%   OUT=matmult(A,B,DIM1) specifies that the matrices are given in DIM1 and
%   DIM1+1. If DIM1 is the last dimensions, matrices are assumed to be
%   given in last two dimensions
%   In this case DIM2=DIM1+1. If DIM1=max(ndims(A),ndims(B)), than
%   DIM2=DIM1 and DIM1=DIM1-1;
%
%   OUT=matmult(A,B,DIM1,DIM2) specifies the two dimensions in which the
%   matrices are given. DIM1 must be smaller than DIM2
%
%   See also: dot, cross, mtimes
%
%   Copyright 2017, Bart Vermeulen

% Check number of inputs
narginchk(2,4); 
nargoutchk(0,1);

% Check properties of input matrices
validateattributes(a,{'numeric'},{},'helpers.matmult','a',1);
validateattributes(b,{'numeric'},{},'helpers.matmult','b',2);

% Number of dimensions involved in computation
nd=max(ndims(a),ndims(b));
if nd == 2
    out = a * b;
    return
end


if nargin > 2                                                              % a dim1 was given
    validateattributes(dim1,{'numeric'},{'integer','scalar','positive',...
        'nonzero'},'helpers.matmult','dim1',3);                               % Check it is ok
    if nargin>3                                                            % a dim2 was given
        validateattributes(dim2,{'numeric'},{'integer','scalar',...
            'positive','nonzero'},'helpers.matmult','dim2',4);                % Check it is ok
        assert(dim1<dim2,'dim1 should be smaller than dim2');              % make sure that dim1<dim1
    else                                                                   % no dim2 was given, so we will compute it
        dim2=dim1+1; 
    end
    assert(dim2<=nd,'dim2 exceeds number of dimensions')                   % make sure dim1 and dim2 do not exceed number of dimensions
else                                                                       % if no dimension is given
    dim1=nd-1;                                                             % dim1 is second to last dimension
    dim2=nd;                                                               % dim2 is last dimension
end

% check inner matrix dimensions agree
assert(size(a,dim2)==size(b,dim1),...
    'Inner dimensions must agree (size(a,dim2)==size(b,dim1))');

% compute output size of matrix product
outsize1=size(a,dim1);
outsize2=size(b,dim2);

% pad size of matrices to both have nd elements (pads with ones for
% trailining singleton dimensions
sizea=[size(a), ones(1,nd-ndims(a))];
sizeb=[size(b), ones(1,nd-ndims(b))];

% Check sizes not involved in multiplication are equal or singleton
sizea_not_mult=sizea;                                                      % Get size of a
sizea_not_mult([dim1,dim2])=[];                                            % Remove multiplication dimensions
sizeb_not_mult=sizeb;                                                      % Get size of b
sizeb_not_mult([dim1,dim2])=[];                                            % Remove multiplication dimensions
assert(all(sizea_not_mult==sizeb_not_mult | sizea_not_mult==1 |...
    sizeb_not_mult==1),...
    ['All dimensions not involved in multiplication must ',...
    'have equal size or must be singleton'])                               % Do the check

% compute size of output
sizeout=max(sizea,sizeb); 
sizeout(dim1)=outsize1;
sizeout(dim2)=outsize2;

% size of matrices that will be multiplied for a and b, i.e. removing the
% non-innner dimensions
sizeamult=sizea;
sizebmult=sizeb;
sizeamult(dim1)=[];
sizebmult(dim2)=[];

% determine permutation order to get right dimensions of matrices to
% multiply
perm_b=1:nd-1;
perm_b(dim1)=dim2-1;
perm_b(dim2-1)=dim1;

% prepare indexing form multiplication
dim1out=shiftdim((1:outsize1)',-dim1+1);                                   % vector counting along the first multiplication dimension
dim2out=shiftdim((1:outsize2)',-dim2+1);                                   % vector counting along the second multiplication dimension

% replication factors to make logical indices match dimensions of matrices
% being referenced
repmat_out=sizeout;                                                        
repmat_out([dim1 dim2])=1;
repmat_a=sizea;
repmat_a(dim1)=1;
repmat_b=sizeb;
repmat_b(dim2)=1;

% compute dimension to sum dot product over
sum_dim=dim2-1; 

% initialize output
out=nan(sizeout,class(a));

% do the computation
for c1=1:outsize1                                                          % loop over dim1
    sel1=dim1out==c1;                                                      % logical indices to select dim1
    for c2=1:outsize2                                                      % loop over dim2
        sel2=dim2out==c2;                                                  % logical indices to select dim2
        out(repmat( sel1 & sel2, repmat_out))=...                          % output value equals
            sum(...                                                        % sum of element product of:
                    reshape(...                                            %        reshaped array
                        a(repmat(sel1,repmat_a)),...                       %            a
                    sizeamult).*...                                        %    and
                permute(...                                                %        permutation to size of a
                    reshape(...                                            %        of the reshaped
                        b(repmat(sel2,repmat_b)),...                       %            b
                    sizebmult),...                                         %
                perm_b),...                                                %        with permutation order perm_b
            sum_dim);                                                      % sum over sum_dim
    end
end