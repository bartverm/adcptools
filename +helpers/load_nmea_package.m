function tf = load_nmea_package()
% Load nmea package
%
%   tf = helpers.load_nmea_package() loads the nmea package if it was not 
%   loaded yet. Returns whether the nmea package is available after the 
%   call to the function.
    tf = nmea_is_loaded();
    if tf
        return
    end
    addpath(fullfile(helpers.adcptools_root, "nmea_submodule"))
    tf = nmea_is_loaded();
end

function tf = nmea_is_loaded()
    tf = exist('nmea.Message','class') == 8;
end