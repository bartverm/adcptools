% This script gives an examble on how to process a vessel mounted adcp
% dataset
restoredefaultpath                                                         % restore default matlab path
clearvars                                                                  % clear all variables
addpath ../../                                                             % point this to adcptools main folder

%% Load adcp data
data=readDeployment('trans','raw_data/');

%% Plot sailed track
[x,y]=utmADCP(data);                                                       % Extract position of ADCP
plot(x,y)                                                                  % Plot track of ADCP
axis equal                                                                 % Set equal axis 

%% Load tid matrix
load tid

%% Process data
msh=procTrans(data,tid,...                                                 % process repeat transects defined in tid
    'CumulateCrossings', true,...                                          % cumlatively average repeat crossings
    'ShipReference','btgps');                                              % use bottom tracking, and if unavailable gps for boat velocity

%% Make simple velocity plot
figure
cs=3;                                                                      % select cross section
rot='cs';                                                                  % select velocity rotation
imagesc(msh(cs).(rot).pars(:,:,end,1))                                     % quickly show longitudinal velocity

%% Make fancy velocity plot (using p structure in output mesh)
figure
cs=3;                                                                      % select cross section
rot='cs';                                                                  % select velocity rotation
qscale=20;                                                                 % set secondary flow quiver scaling


patch(msh(cs).p.N,...                                                      % n-coordinate of cell edges
      msh(cs).p.Z(:,:,end),...                                             % z-coordinate of cell edges (end: selects average over all repeats)
      msh(cs).(rot).pars(msh(cs).p.progfgood_vec(:,end,1)))                % longitudinal flow (used progfgood to select only values inside mesh
shading flat                                                               % remove edges of cells
set(gca,'clim',[-0.3 0.6])                                                 % set color limits of figure
colormap(flipud([103,0,31;                                                 % set a colormap (based on colorbrewer)
178,24,43;  
214,96,77;
244,165,130;
253,219,199;
247,247,247;
247,247,247;
209,229,240;
146,197,222])/255)
cb=colorbar;                                                               % add a colorbar
hold on                                                                    % add next plots
plot(msh(cs).p.nbed,msh(cs).p.zbed,'k','linewidth',2)                      % draw bed
plot(msh(cs).p.nbed([1 end]),[0 0],'b','linewidth',2)                      % draw water surface
quiver(...                                                                 % draw secondary flow vectors
       msh(cs).N,...                                                       % n coordinate of cell centers
       msh(cs).Z(:,:,end),...                                              % z coordinate of cell centers
       squeeze(msh(cs).(rot).pars(:,:,end,2)*qscale),...                   % n-velocity
       squeeze(msh(cs).(rot).pars(:,:,end,3)*qscale),...                   % z-velocity
       'color','k','autoscale','off')                                      % black quivers, manually scaled

% add labels and quiver scale to plot
xlabel('n coordinate (m)')
ylabel('z coordinate (m)')
ylabel(cb,'Longitudinal velocity (m/s)')
quiver(50, -35,1*qscale,0,'k','autoscale', 'off');
quiver(50, -35, 0,0.2*qscale,'k','autoscale', 'off');
text(45,-33,'0.2 m/s','HorizontalAlignment','right')
text(60,-36,'1 m/s','HorizontalAlignment','center',...
                    'VerticalAlignment','top')

