function out=matmult(A,B,dims)
% Matrix multiplication

% dims must be 2 element vector with integers
assert(isvector(dims) & numel(dims)==2 & all(floor(dims)==dims));

% Inner matrix dimensions must agree
assert(size(A,dims(2))==size(B,dims(1)));

sizeA=size(A);
sizeB=size(B);
if numel(sizeA) < numel(sizeB)
    sizeA=[sizeA ones(1,numel(sizeB)-numel(sizeA))];
elseif numel(sizeA) > numel(sizeB)
    sizeB=[sizeB ones(1,numel(sizeA)-numel(sizeB))];
end

outsize=max([sizeA;sizeB],[],1);
m=sizeA(dims(1));
n=sizeB(dims(2));
% n=size(A,dims(2));
outsize(dims)=[m n];
out=nan(outsize);
totdims=numel(outsize);

evstr=':,';
Astr=[repmat(evstr,1,max(dims(1)-1,0)),'cr,',repmat(evstr,1,totdims-dims(2)+1)];
Bstr=[repmat(evstr,1,dims(1)),'cc,',repmat(evstr,1,totdims-dims(2))];
Outstr=[repmat(evstr,1,max(dims(1)-1,0)),'cr,cc,',repmat(evstr,1,totdims-dims(2))];
Outstr(end)=[];Outstr=['out(',Outstr,')=sum(bsxfun(@times,Acur,Bcur),dims(2));'];
Astr(end)=[];Astr=['A(',Astr,')'];
Bstr(end)=[];Bstr=['B(',Bstr,')'];
perm=1:totdims;
perm(dims)=perm([dims(2) dims(1)]);

for cr=1:m
    for cc=1:n
        Acur=eval(Astr); %#ok<NASGU>
        Bcur=permute(eval(Bstr),perm); %#ok<NASGU>
        eval(Outstr);
    end
end
