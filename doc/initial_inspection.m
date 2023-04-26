%% Initial inspection of ADCP data
% The classes <matlab:doc('ADCP') ADCP> and <matlab:doc('VMADCP') VMADCP>
% implement a number of methods that help the user to inspect the data. 
% These functions are mostly plot functions
%

%% Reading the data
% We start by reading the data:

%%
% raw_dat = rdi.readDeployment('trans',...
%     fullfile(helpers.adcptools_root, 'doc','sample_data',...
%     'rdi_muara_muntai_bend'));
%
%%
% Next we construct a VMADCP object:

%%
v = rdi.VMADCP(raw_dat);

%% Instrument orientation
% Insepct the tilts of the instrument during deployment:

%%
v.plot_orientations

%%
% The orientations are within reasonable limits for ship borne deployments

%% Ship track
% We can plot the ship track as follows:

%%
figure
v.plot_track

%% Detected bed elevation
% The bed elevations as detected by each of the acoustic beams and
% corrected for tilts of the instrument can be displayed with:

%%
figure
v.plot_bed_position

%%
% This plot displays the elevation of the bed. For this computation to
% properly work the elevation of the ADCP needs to be known. These data
% were collected in a river with negligible elevation changes. The way the
% vertical position is determined, is defined in the 
% <matlab:doc('ADCP.vertical_position_provider') vertical_position_provider>
% property of the ADCP object.

%%
v.vertical_position_provider

%%
% The vertical position is provided by an
% <matlab:doc('ADCPVerticalPositionFromWaterLevel') ADCPVerticalPositionFromWaterLevel>
% object. This object has a
% <matlab:doc('ADCPVerticalPositionFromWaterLevel.water_level') water_level>
% property that defines the waterlevel. This is currently set to 
% <matlab:doc('ConstantWaterLevel') ConstantWaterLevel> object. Furthermore
% the 
% <matlab:doc('ADCPVerticalPositionFromWaterLevel.depth_transducer') depth_transducer> 
% property defines how deep the transducer is inserted in the water. For a
% depth of 40 cm we would set this property as:

%%
v.vertical_position_provider.depth_transducer = 0.4;

%%
% We can take a further look at the 
% <matlab:doc('ADCPVerticalPositionFromWaterLevel.water_level') water_level>
% property

%%
v.vertical_position_provider.water_level

%% 
% The <matlab:doc('ConstantWaterLevel') ConstantWaterLevel> object has the
% property <matlab:doc('ConstantWaterLevel.level') level> set to 0. This
% means the water level is at an elevation of  0 m. You may want to change
% this to have all position computations with respect to a given datum.
% Suppose the water level was at 60 m above sea level we can set:

%%
v.vertical_position_provider.water_level.level = 60;

%%
% We will see the changes reflected in the bed elevation computation:

%%
figure
v.plot_bed_position;

%% Depth averaged velocity
% The depth averaged velocity is plotted as:

%%
figure
v.plot_track_velocity


%% Velocity data
% The velocity data can be inspected as follows:

%%
figure
v.plot_velocity
set(gca,'xlim',[1500 2500]) % zooming in on portion of the data

%%
% Note that by default the velocity is shown in earth coordinate system,
% i.e. Vx is east velocity, Vy is north velocity and Vz is upward velocity.








