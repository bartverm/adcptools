function varargout = readDeployment(varargin)
    warning('readDeployment:obsolete',...
        ['This function will be removed in the future, ',...
        'use rdi.readDeployment instead.'])
    varargout = cell(nargout);
    [varargout{:}] = rdi.readDeployment(varargin{:});
end