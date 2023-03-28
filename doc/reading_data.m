%% Reading data
% How to read data
%
%% Introduction
% The starting point to analyze ADCP data is to load the ADCP data files.
% The way to read the ADCP data depends on the manufacturer and kind of
% dataset you want to read in. ADCPtools currently supports the following
% manufacturers
%
% * RDI
% * Nortek
% * Sontek
%
% In order to further process the data, it is necessary to construct an
% <matlab:doc('ADCP') ADCP> object for self-contained data and a <matlab:doc('VMADCP') VMADCP> objects for
% vessel-mounted data. These classes contain several functions that allow
% to process and visualize the data.

%% Reading RDI data
% ADCP tools always requires a raw, binary PD0 file to be read. For RDI
% data, the reading of data and construction of <matlab:doc('ADCP') ADCP> or <matlab:doc('VMADCP') VMADCP> objects is
% done in two separate steps. First the data are read in a structure, and
% subsequently, the <matlab:doc('ADCP') ADCP> or VMADCP objects are constructed from these
% structures.
% To read in a self-contained/moored dataset you can use the rdi.readADCP
% function:
%
%   dat = rdi.readADCP('filename.PD0');
% 
% This will read the PD0 data and return a data structure containing all
% data in the PD0 file. After reading the raw data, an ADCP object can be
% constructed:
% 
%   Dataset = rdi.ADCP(dat);
%
% This constructs an <matlab:doc('rdi.ADCP') rdi.ADCP> object, which is a subclass of <matlab:doc('ADCP') ADCP>.
%
% When reading data from a vessel mounted deployment, next to the raw PD0
% files you will also have external data such as GPS, echo-sounder data or
% external heading. These data are usually organized in set of files
% starting with a common part called the deployment name followed by a
% number which identifies the repeat transect. To read an entire deployment
% you can use the rdi.readDeployment function:
%
%   dat = rdi.readDeployment('deploymentName','path/')
%
% This will read all data files starting with 'deploymentName'. The
% optional argument 'path' allows to specify the location of the data. The
% next step is to create a <matlab:doc('VMADCP') VMADCP> object:
%
%   VMDataset = rdi.VMADCP(dat);

%% Reading Nortek data
% The ADCPtools currently supports reading AD2CP files for self-contained
% deployments and SigVM file for vessel mounted deployments. For
% self-contained deployments data can be read as:
%
%   Dataset = nortek.ADCP('filename.ad2cp');
%
% This reads all data from the ad2cp file and constructs the <matlab:doc('ADCP') ADCP> object
% for further processing and visualization of the data.
%
% For a vessel mounted dataset, data is read as follows:
%
%   VMDataset = nortek.VMADCP('filename.SigVM');
%
% This reads the ad2cp data and the gnss data in the SigVM file and
% constructs the <matlab:doc('VMADCP') VMADCP> object

%% Reading Sontek data
% ADCPtools can construct <matlab:doc('ADCP') ADCP> and <matlab:doc('VMADCP') VMADCP> objects based on the matlab
% exports which can be generated with the RiverSurveyor software from
% Sontek. Once the .mat files have been generated they can be read as
% follows:
%
%   Dataset = sontek.VMADCP('filename.mat');
%
% This constructs the <matlab:doc('VMADCP') VMADCP> object for further processing.
