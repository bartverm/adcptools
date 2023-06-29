%% Solving repeat transect data on the mesh

%%
% The last step in the processing of the repeat-transect data is the
% solution of the velocity field. This step is handled by the 
% <matlab:doc('VelocitySolver') VelocitySolver class>. This is a generic
% class that cannot directly be used. Two implementations exist:
%%
%
% * <matlab:doc('TimeBasedVelocitySolver') TimeBasedVelocitySolver>: 
% Implements the classic
% way of processing ADCP data in which beam velocities are combined to
% obtain a Cartesian velocity based on the time at which they were measured
% * <matlab:doc('LocationBasedVelocitySolver) LocationBasedVelocitySolver>:
% is a solver
% that combines the beam velocities based on where they were measured in
% space. This reduces the spatial homegeneity assumption
%
%%
% In case you doubt use the first one. The second solver is recommended
% when measuring in highly sheared flows, and when positioning is of good
% quality. We will use the TimeBasedVelocitySolver:
solver = TimeBasedVelocitySolver(mmbend, mesh, xs, ef, B);

%%
% This again produces seven solvers, one for each cross-section. The last
% step is to solve for the velocity:
vel = solver.get_velocity()

%%
% This returned a 1x7 cell each containing a matrix of Nx3, with N being
% the number of cells in the mesh, and 3 the velocity components. The
% velocity is solved in East-North-Up components. We may want to rotate the
% velocity to the direction of the cross-sections:
vel_sn = solver.rotate_to_xs(vel);

%% 
% We can plot the velocity as follows:
figure
mesh.plot(vel_sn)

%%
% Or in 3D:
figure
mesh.plot3(vel_sn)
