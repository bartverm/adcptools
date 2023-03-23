classdef ShipBubbleFilter < Filter
   properties
       max_search_depth(1,1) double {mustBeNonnegative} = 3;    % m from surface
       max_bs_deviation(1,1) double {mustBeNonnegative} = 10;    % mark cells as bad when backscatter value deviates from median horizontal backscatter with max_bs_deviation
   end
   methods
        function obj=ShipBubbleFilter()
            obj.description='Ship bubble filter';
        end
    end
   methods (Access=protected)
       function bad=bad_int(obj,adcp)
           validateattributes(adcp,{'VMADCP'},{'scalar'})
           if isscalar(obj)
               bad = zeros(size(adcp.echo));
               BS = adcp.backscatter(Filter);
               bs_median = median(BS,2,'omitnan');   % median backscatter in horizontal direction


               binsize = adcp.cellsize(1);
               n_topbins = round(obj.max_search_depth/binsize);
               maxdev = obj.max_bs_deviation;

               % Start marking bad cells near the surface
               for nbin = 1:n_topbins
                   for b = 1:4         % Loop over beams
                       remidx = BS(nbin,:,b) - bs_median(nbin,1,b) > maxdev;
                       bad(nbin,remidx,b) = 1;
                   end
               end

               % For lower layers, only filter if cell on top/neighbouring
               % on top has also been removed (this part is very slow)
               for nbin = n_topbins+1:size(BS,1)
                   for b = 1:4         % Loop over beams
                       [~, badcolls] = find(bad(nbin-1,:,b) == 1);  % Bad cell indices
                       for n = 1:length(badcolls)
                           c = badcolls(n);
                           if BS(nbin,c,b) - bs_median(nbin,1,b) > maxdev && ((bad(nbin-1,c,b)==1) || (bad(nbin-1,c-1,b)==1) || (bad(nbin-1,c+1,b)==1))
                               bad(nbin,c,b) = 1;
                           end
                       end
                   end
               end
            bad = logical(bad);
           else
               bad=bad_int@Filter(obj,adcp);
           end
       end
   end
end
