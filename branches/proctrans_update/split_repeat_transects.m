function ef_out=split_repeat_transects(v,xs,varargin)
% splits data in a cross-section to different repeat transects
%
%   ef_out = split_repeat_transects(v,xs) returns a vector of
%   ensemblefilter object. Each object indicates the ensembles in the
%   VMADCP object v that belong a specific repeat transect on the 
%   cross-section xs. The algorithm detects peaks in the n-coordinate of 
%   the track, when the peak is spaced at least the scale of the 
%   cross-section (xs.scale) from the previous peak.
%
%   split_repeat_transects(..., scale) optionally specify scale for peak
%   detection algorithm.
%
%   split_repeat_transect(..., ef) optionally specify an ensemble filter
%   object to exclude part of the data
%
%   see also: cross_section_selector, helpers.peakdet, VMADCP,
%   CrossSection, EnsembleFilter

scale=xs.scale;
ef=EnsembleFilter(v);
for ca=1:numel(varargin)
    cur_arg=varargin{ca};
    if isa(cur_arg,'double') && isscalar(cur_arg) && isfinite(cur_arg)
        scale=cur_arg;
    elseif isa(cur_arg,'EnsembleFilter') && isscalar(cur_arg)
        ef=cur_arg;
    else
        warning(['Unhandled input number ', num2str(ca+2), 'of type: ', class(cur_arg)])
    end
end

fgood=~ef.bad_ensembles; % get which ensembles to include
hpos=v.horizontal_position; % get horizontal position of ADCP
[~,n]=xs.xy2sn(hpos(1,fgood),hpos(2,fgood)); % transform to direction along cross-section
[maxidx,minidx]=helpers.peakdet(n,scale); % detect minima and maxima in n-coordinate with scale of the cross-section
idx=sort([maxidx(:,1); minidx(:,1)]); % lump and sort minima and maxima

nrp=size(idx,1)-1; % number of repeat transects
ef_out(nrp,1)=EnsembleFilter; % initialize output ensemble filter array
fgood=find(fgood); % get indices of ensembles included in repeat transect detection (to transform peakdet indices to adcp indices)
for cc=1:nrp % loop over all repeat transects
    cfilt=true(1,v.nensembles); % initialize filter to filter out eveything
    cfilt(fgood(idx(cc)):fgood(idx(cc+1)))=false; % set current repeat to false
    cfilt=cfilt | ef.bad_ensembles; % make sure no ensembles get included that where filtered out by given ensemble filter
    ef_out(cc)=EnsembleFilter(v,cfilt); % create output ensemble filter object.
end