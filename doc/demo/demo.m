% This script gives an examble on how to process a vessel mounted adcp
% dataset
restoredefaultpath                                                         % restore default matlab path
clearvars                                                                  % clear all variables
addpath ../../                                                             % point this to adcptools main folder
addpath ~/src/loess/ % For bathymetry interpolation, download at https://github.com/bartverm/loess

%% Load adcp data
data=readDeployment('trans','raw_data/');

%% Create VMADCP object
V = rdi.VMADCP(data);

%% Make some plots (for first inspection of data)
V.plot_all


%% Define cross sections
load xsdefinitions % load predefined cross-sections
% or
% [ef, xs] = cross_section_selector(V); % interactive tool to aid definition of cross sections


%% Create the meshes for velocity solutions
for cs = 1 : numel(ef)
    Bathy(cs) = BathymetryScatteredPoints(V,ef(cs)); % This can also be done for each cross-section separately
    Bathy(cs).interpolator.span=0.005; % lower the degree of smoothing (neighborhood is now 0.5 % of all depth points)
    MMaker(cs) = SigmaZetaMeshFromVMADCP(V,...
        Bathy(cs),...
        xs(cs),...
        ef(cs)); % object to generate meshes
    MMaker(cs).deltan = 5; % lateral target resolution
    MMaker(cs).deltaz = 1; % vertical target resolution
    Mesh(cs) = MMaker(cs).get_mesh();
    Vsolver(cs) = LocationBasedVelocitySolver(V,... % Use TimeBasedVelocitySolver to solve in the conventional way, i.e. using homogeneity assumption
        Bathy(cs),...
        xs(cs),...
        ef(cs),...
        Mesh(cs));
    [vel(cs), cov_vel(cs)] = Vsolver(cs).get_velocity(); % solve the velocity and also get covariances in velocity
    [vel_sn(cs), cov_vel_sn(cs)] = Vsolver(cs).rotate_to_xs(vel(cs),cov_vel(cs)); % Rotate velocity and covariances to cross-section direction
end

%% Make some plots
for cs = 1:numel(ef)
    % Bathymetry and mesh
    figure
    title('Bathymetry and mesh');
    Bathy(cs).plot;
    hc = colorbar;
    ylabel(hc, 'bed elevation (m)')
    hold on
    xs(cs).plot;
    Mesh(cs).plot3;
    xlabel('UTM x (m)')
    ylabel('UTM y (m)')
    zlabel('Elevation (m)')

    figure
    title('Velocity')
    Mesh(cs).plot(vel_sn{cs}(:,1))
    hc = colorbar;
    ylabel(hc, 'Downstream velocity (m/s)')
end


