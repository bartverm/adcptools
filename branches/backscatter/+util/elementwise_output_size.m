function siz=elementwise_output_size(varargin)
% Returns output size of elementwise operation on input arrays
%
% SIZ=elementwise_output_size(A1, A2, A3, ...) outputs the size of an 
%   output array resulting from an elementwise computation with the input \
%   arrays. Gives an error if dimensions do not match

nd=max(cellfun(@ndims,varargin));
sizes=cellfun(@(x) [size(x) max(0,ones(1,nd-ndims(x)))],varargin,'UniformOutput',false);
sizes=vertcat(sizes{:});
siz=max(sizes,[],1);
assert( all(all(sizes==siz | sizes==1)),'elementwise_output_size:incompatible','Arrays are incompatible for elementwise computation')