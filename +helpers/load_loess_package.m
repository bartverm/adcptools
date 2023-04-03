function tf = load_loess_package()
% Load loess package
%
%   tf = helpers.load_loess_package() loads the nmea package if it was not 
%   loaded yet. Returns whether the nmea package is available after the 
%   call to the function.

    % Case loess exists and runs
    tf = loess_exists && loess_runs;
    if tf
        return
    end

    adcptools_path = fileparts(which('ADCP'));
    loess_path = fullfile(adcptools_path, "loess_submodule");


    % If loess does not exist, add path
    if ~loess_exists
        addpath(loess_path);
    end

    % check if loess exists and if it runs
    tf = loess_exists && loess_runs;
end


function tf = loess_exists()
    tf = exist(['loess.', mexext],'file') == 3;
end

function tf = loess_runs()
    try
        loess
    catch err
        tf = strcmp(err.identifier,'loess:nargin');
    end
end
