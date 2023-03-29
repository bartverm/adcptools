function tf = load_loess_package()
% Load loess package
%
%   tf = helpers.load_loess_package() loads the nmea package if it was not 
%   loaded yet. Returns whether the nmea package is available after the 
%   call to the function.
    tf = loess_is_loaded();
    if tf
        return
    end
    adcptools_path = fileparts(which('ADCP'));
    loess_path = fullfile(adcptools_path, "loess_submodule");
    addpath(loess_path)
    tf = loess_is_loaded();
    if ~tf && exist("compile_loess.m",'file') == 2
        current_path = pwd;
        cd(loess_path)
        try
            compile_loess;
        catch
            warning('Tried to compile loess but failed.')
        end
        cd(current_path);
        tf = loess_is_loaded();
    end

end

function tf = loess_is_loaded()
    tf = exist(['loess.', mexext],'file') == 3;
end