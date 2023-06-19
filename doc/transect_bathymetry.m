%% Creation of a bathymetric model

%%
% Once the VMADCP data were selected and cross-sections are defined, it is
% necessary to define a bathymetric model, which will be used in the rest
% of the processing. The most obvious way, is to construct the bathymetry
% model from the VMADCP bottom tracking data. In some cases, however, these
% data maybe of poor quality and a bathymetric model might be constructed
% from another source.

%% The Bathymetry class
% The requirements for a bathymetric model are defined in the 
% <matlab:doc('Bathymetry') Bathymetry> class. In essence a 
% <matlab:doc('Bathymetry') Bathymetry> object must be able to provide a
% bed elevation for any given (x,y) position.

%% Bathymetry from VMADCP data
% The Bottom tracking of ADCPs produces a scattered cloud of bed
% detection. The 
% <matlab:doc('BathymetryScatteredPoints') BathymetryScatteredPoints> class
% defines a bathymetric model from a given clous of scattered points. The
% class used an interpolator to smooth and interpolate the input cloud.
% The <matlab:doc('BathymetryScatteredPoints') BathymetryScatteredPoints>
% can be constructed from a <matlab:doc('VMADCP') VMADCP> object. In that
% case the class obtains the bed positions from the ADCP bottom tracking
% and uses those as source for the cloud of bed positions.
% In the case of the Muara Muntai bend dataset we have seven
% cross-sections. For each we wish to construct a bathymetry. We can do
% this by constructing the 
% <matlab:doc('BathymetryScatteredPoints') BathymetryScatteredPoints> by
% passing it the <matlab:doc('VMADCP') VMADCP> object and the 
% <matlab:doc('EnsembleFilter') EnsembleFilter> objects:
B = BathymetryScatteredPoints(mmbend, ef);


%% Creating a custom Bathymetry
