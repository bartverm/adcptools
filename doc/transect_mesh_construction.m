%% Mesh creation for repeat transect data processing

%%
% The next step in the processing of ADCP data is the creation of a mesh on
% which we will eventually solve the velocity. The mesh that we will create
% is a so-called sigma-zeta mesh. It is a kind of hybrid between a zeta
% mesh, in which vertical cells all have equal size and a sigma mesh where
% the number of cells is everywhere the same. The mesh is implementend in
% the <matlab:doc('SigmaZetaMesh') SigmaZetaMesh> class, which can be
% constructed with a 
% <matlab:doc('SigmaZetaMeshGenerator') SigmaZetaMeshGenerator>. In this 
% case we will use the 
% <matlab:doc('SigmaZetaMeshFromVMADCP') SigmaZetaMeshFromVMADCP> class.
% The class needs objects of the following classes
%%
% 
% * <matlab:doc('VMADCP') VMADCP>: to provide details about the flow
% measurements
% * <matlab:doc('Bathymetry') Bathymetry>: to know the bed position
% * <matlab:doc('EnsembleFilter') EnsembleFilter>: to know which data to
% include
% * <matlab:doc('XSection') Xsection>: for the definition of a cross
% section.

%% Mesh constructor creation
% We first construct the mesh generator object:
mesh_gen = SigmaZetaMeshFromVMADCP(mmbend, ef, xs, B)

%%
% We get one generator for each cross-section. Some properties of the mesh
% generator are useful to consider:
%%
%
% * _deltan_: lateral resolution of the mesh
% * _deltaz_: vertical resolution of the mesh
% * _water_level_: defines water level
% * _time_: time for which mesh is to be generated
%
%%
% The first two properties are quite straightforward, but the last two
% might need some explanation. The water_level property is important
% particularly when dealing with situations in which water level variations
% are significant during the measurements. In this case, also the time for
% which the mesh should be generated is important, since it sets the water
% level for which it will be generated.

%% Generate the meshes
% We can now generat the mesh
mesh = mesh_gen.get_mesh();

%% 
% Seven meshes were generated, one for each cross-section. We can visualize
% the meshes.
figure
mesh.plot

%% 
% we can also visualize the meshes in 3D:
figure
B.plot
hold on
mesh.plot3
