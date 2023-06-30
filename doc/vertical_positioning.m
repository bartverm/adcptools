%% Vertical positioning
% An important step when processing ADCP data is to properly set the
% vertical positioning of the data. For moored ADCP deployments the
% vertical position is usually set to a fixed level. For vessel mounted
% deployments the vertical position of the ADCP follows the water surface.
% This water surface might be (quasi) constant in a river, but might vary
% significantly in a tidally affected area. 

%% Moored ADCP: fixed vertical position


%% Vessel mounted ADCP: constant water surface
% If we first take a look at the velocity measured by the ADCP:

%%
mmbend.plot_velocity
set(gca, 'xlim', mmbend.time([1500 2500])) % zooming in on portion of the data

%%
% We can see that by default the water surface is at 0 m elevation, and
% this position is constant during the entire measurement. These data
% were collected in a river with negligible elevation changes. The way the
% vertical position is determined, is defined in the 
% <matlab:doc('ADCP.vertical_position_provider') vertical_position_provider>
% property of the ADCP object.

%%
mmbend.vertical_position_provider

%%
% The vertical position is provided by an
% <matlab:doc('ADCPVerticalPositionFromWaterLevel') ADCPVerticalPositionFromWaterLevel>
% object. This object uses the 
% <matlab:doc('ADCP.water_level') water_level> and
% <matlab:doc('VMADCP.depth_transducer') depth_transducer> properties of the 
% <matlab:doc('ADCP') ADCP> object to compute the vertical position of the ADCP.
% We set the depth of the transducer to 0.4 m:
%%
mmbend.depth_transducer = 0.4;

%%
% The <matlab:doc('ADCP.water_level') water_level> property of the 
% <matlab:doc('ADCP') ADCP> object is controlled by the 
% <matlab:doc('ADCP.water_level_object') water_level_object> property

%%
mmbend.water_level_object

%% 
% The <matlab:doc('ConstantWaterLevel') ConstantWaterLevel> object has the
% property <matlab:doc('ConstantWaterLevel.level') level> set to 0. This
% means the water level is at an elevation of  0 m. You may want to change
% this to have all position computations with respect to a given datum.
% Suppose the water level was at 60 m above sea level we can set:

%%
mmbend.water_level_object.level = 60;

%%
% We will see the changes reflected in all position computations of the
% ADCP:

%%
figure
mmbend.plot_velocity;
set(gca, 'xlim', mmbend.time([1500 2500])) % zooming in on portion of the data

%%
% we set the vertical position back to 0:
mmbend.water_level_object.level = 0;

%% Vessel mounted ADCP: varying water surface

%% Customize vertical positioning
