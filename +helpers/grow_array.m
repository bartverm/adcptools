function varargout=grow_array(varargin)
% Grow arrays to have matching dimensions
%
%   out=grow_array(a1, a2, ..., an) Grow arrays to have matching dimensions
%
%   out=fillcat(...,'Name',value) specify additional options with
%   name-value pairs:
%
%   'FillVal' 
%   specify the value to grow arrays with. For handle arrays, the object is
%   grown with unique new objects. If 'FillVal' is specified, it will be
%   filled with copies of 'FillVal'
%
%   'SkipDims'
%   specify dimensions that should not be grown. Default is empty If both
%   SkipDims and OnlyDims are given, SkipDims is ignored.
%
%   'OnlyDims'
%   only grow arrays along given dimensions. Default is empty. If both
%   SkipDims and OnlyDims are given, SkipDims is ignored.
%
%   'CollectOutput'
%   Collect all grown arrays in a cell. Default is false


% *** Parse input ***
if nargin < 1
    return;
end

if nargin < 2
    varargout{1}=varargin{1};
    return
end

P=inputParser;
P.addParameter('FillVal',make_fillval(varargin{1}),@(x) isscalar(x) && isequal(class(x),class(varargin{1})));
P.addParameter('SkipDims',[],@(x) isnumeric(x) && all(x>0) && all(mod(x,1)==0));
P.addParameter('OnlyDims',[],@(x) isnumeric(x) && all(x>0) && all(mod(x,1)==0));
P.addParameter('CollectOutput',false,@(x) islogical(x) && isscalar(x))

% Compute number of input variables
n_var=1; % we start with first
while n_var < nargin &&... % if we still have unprocessed inputs
      isequal(class(varargin{n_var}), class(varargin{n_var+1})) &&... % and their class matches the class of the previous input
      ~(ischar(varargin{n_var+1}) && any(strcmp(varargin{n_var+1},P.Parameters))) % And it is not an parameter name - Note that this excludes the possibility to concatenate parameter names
    n_var=n_var+1; % increase the number of inputs
end
if n_var==1 % If only one variable of same class was found
    error('At least the first two inputs must be of the same class') % generate an error
end

P.parse(varargin{n_var+1:nargin})
fillval=P.Results.FillVal;
skipdims=P.Results.SkipDims;
onlydims=P.Results.OnlyDims;
collect_output=P.Results.CollectOutput;

fillsiz=size_equal(varargin(1:n_var));
target_siz=max(fillsiz,[],1);
n_grow=max(target_siz-fillsiz,0);
dims=find(any(n_grow~=0,1));
[~, skipidx] = intersect(dims,skipdims);
dims(skipidx)=[];
if ~isempty(onlydims)
    dims = reshape(onlydims,1,[]);
end
for cd=dims
     vars=reshape(find(n_grow(:,cd)>0),1,[]);
     for cvar=vars
         fillsiz(cvar,cd)=n_grow(cvar,cd);
         if isa(varargin{1},'handle') && ~isempty(strcmp('FillVal',P.UsingDefaults))
             tmp=make_handle_fill(fillsiz(cvar,:),class(varargin{1}));
         else
             tmp=repmat(fillval,fillsiz(cvar,:));
         end
         varargin{cvar}=cat(cd,varargin{cvar},tmp);
     end
     fillsiz=size_equal(varargin(1:n_var));     
end

if collect_output
    varargout{1}=varargin(1:n_var);
else
    varargout=varargin(1:n_var);
end

end

function out=make_handle_fill(siz,class_name) %#ok<STOUT>
    evstr=num2str(reshape(siz,[],1));
    evstr=[evstr, repmat(',',size(evstr))];
    evstr=reshape(evstr',1,[]);
    evstr(end)=[];
    if any(siz==0)
        eval(['out=',class_name,'.empty(',siz,');'])
    else
        eval(['out(',evstr,')=',class_name,';'])
    end 
end
function out=make_fillval(a)
    if isfloat(a)
        out = cast(NaN,'like', a);
        return
    end
    if isinteger(a)
        out = intmin(class(a));
        return;
    end
    if ischar(a)
        out = ' ';
        return
    end
    try
        out=cast(0,'like',a);
    catch
        out=eval(class(a));
    end
end 
function siz=size_equal(in)
nd=max(cellfun(@ndims,in));
siz=cellfun(@(x) [size(x) ones(1,nd-ndims(x))],in,'UniformOutput',false);
siz=vertcat(siz{:});
end