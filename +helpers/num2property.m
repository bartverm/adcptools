function num2property(obj, property, var, dim)
% Assign array values to property of an object array
%
%   num2property(OBJ, PROPERTY, VAR) assigns the values in VAR to the
%   PROPERTY of the object in the array OBJ. This is done along the
%   dimension with the same size as number of objects in OBJ. 
%
%   num2property(OBJ, PROPERTY, VAR, DIM) also specify the dimension along
%   which to split the varibale

if isscalar(var)
    var=repmat(var,size(obj));
end
if nargin < 4
    nobj=numel(obj);
    siz_var=size(var);
    dim_obj=find(siz_var==nobj);
    dim=find(siz_var~=nobj);
    assert(numel(dim_obj)==1 && numel(dim)==1,...
  'cannot find dimension in var with same number of elements as the object array')
end
tmp=num2cell(var,dim);
[obj.(property)]=deal(tmp{:});


end