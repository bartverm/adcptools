function [maxtab, mintab]=peakdet(v, delta, x)
%PEAKDET Detect peaks in a vector
%        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
%        maxima and minima ("peaks") in the vector V.
%        MAXTAB and MINTAB consists of two columns. Column 1
%        contains indices in V, and column 2 the found values.
%      
%        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
%        in MAXTAB and MINTAB are replaced with the corresponding
%        X-values.
%
%        A point is considered a maximum peak if it has the maximal
%        value, and was preceded (to the left) by a value lower by
%        DELTA.

% Eli Billauer, 3.4.05 (Explicitly not copyrighted).
% This function is released to the public domain; Any use is allowed.

%maxtab = [];
%mintab = [];

v = v(:); % Just in case this wasn't a proper vector

if nargin < 3
  x = (1:length(v))';
else 
  x = x(:);
  if length(v)~= length(x)
    error('Input vectors v and x must have same length');
  end
end
  
if (length(delta(:)))>1
  error('Input argument DELTA must be a scalar');
end

if delta <= 0
  error('Input argument DELTA must be positive');
end

mn = Inf; mx = -Inf;
mnpos = NaN; mxpos = NaN;

lookformax = 1;
maxtab = nan(length(v),2);
mintab = nan(length(v),2);

for i=1:length(v)
  this = v(i);              % Set this to current value
  if this > mx, mx = this; mxpos = x(i); end  % If current value is higher than mx, change mx to current value and store current position
  if this < mn, mn = this; mnpos = x(i); end  % If current value is smaller than minimum, change mn to current value and store current position
  
  if lookformax
    if this < mx-delta
      %maxtab = [maxtab ; mxpos mx];
      maxtab(i,:) = [mxpos mx];
      mn = this; mnpos = x(i);
      lookformax = 0;
    end
  else
    if this > mn+delta
      %mintab = [mintab ; mnpos mn];
      mintab(i,:) = [mnpos mn];
      mx = this; mxpos = x(i);
      lookformax = 1;
    end
  end
end
maxtab = maxtab(~isnan(maxtab(:,1)),:);
mintab = mintab(~isnan(mintab(:,1)),:);