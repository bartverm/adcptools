function varargout = readADCP(varargin)
    warning('readDeployment:obsolete',...
        ['This function will be removed in the future, ',...
        'use rdi.readADCP instead.'])
    varargout = cell(nargout);
    [varargout{:}] = rdi.readADCP(varargin{:});
end