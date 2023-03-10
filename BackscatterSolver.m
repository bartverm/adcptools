classdef BackscatterSolver < ADCPDataSolver
    methods (Access = protected)
        function [vpos, vdat, xform, time, wl] = get_solver_input(obj)
            [vpos, ~, ~, time, wl] = get_solver_input@ADCPDataSolver(obj);             % call superclass method to get positions and time of adcp data
          
            % get backscatter data
            vdat = obj.adcp.backscatter;

            % get transformation matrix
            xform = ones([1, size(vdat,2), size(vdat,3)]);

            % filter and vectorize
            [vdat, xform] = obj.filter_and_vectorize(vdat, xform);
        end
    end
end