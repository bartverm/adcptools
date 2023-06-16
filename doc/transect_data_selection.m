%% Selection of data for repeat transect processing
% When processing repeat transect ADCP data, a selection needs to be made
% regarding which data belongs to a certain cross-section, we also may want
% to exclude parts of the data from the processing. Data selection is done
% through <matlab:doc('EnsembleFilter') EnsembleFilter> objects. An 
% <matlab:doc('EnsembleFilter') EnsembleFilter> object. We can create an
% object:

%%
ef = EnsembleFilter(mmbend)

%% 
% The <matlab:doc('EnsembleFilter') EnsembleFilter> object has the property 
% <matlab:doc('EnsembleFilter.bad_ensembles') bad_ensembles> that marks
% which ensembles are not to be included in the processing. There are
% different ways to generate an 
% <matlab:doc('EnsembleFilter') EnsembleFilter> object.

%% File based selection
% Sometimes during data collection separate files are made for different
% cross-sections or repeat transects. In this case it is possible to select
% data based on the file. The example dataset was generated from two files:

%%
unique(mmbend.fileid)

%%
% We can now create an <matlab:doc('EnsembleFilter') EnsembleFilter> object
% that excludes data that was not in the first file:
ef = EnsembleFilter(mmbend, mmbend.fileid ~= 1);

%%
% we can plot the result:

%%
figure
plot(ef,mmbend)

%% 
% We can see that the last part of the ensembles is now marked as bad and
% will be excluded from the processing. In this case, the selection based
% on the input file does not make sense and we want to use another method.

%% Location based selection
% in order to select the data based on the location we can use the helper
% function <matlab:doc('cross_section_selector') cross_section_selector>. 
% This will present the use with a plot of the track and will allow to draw
% a polygon around the data that will be included in the processing.
% Multiple cross-sections can be selected with this tool



