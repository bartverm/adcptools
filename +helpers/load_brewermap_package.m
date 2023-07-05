function tf = load_brewermap_package()
% Load brewermap package
%
%   tf = helpers.load_brewermap_package() loads the brewermap package if it
%   was not loaded yet. Returns whether the brewermap package is available 
%   after the call to the function.
    tf = brew_is_loaded();
    if tf
        return
    end
    addpath(fullfile(helpers.adcptools_root, "brewermap_submodule"))
    tf = brew_is_loaded();
end

function tf = brew_is_loaded()
    tf = exist('brewermap','file') == 2;
end