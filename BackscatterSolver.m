classdef BackscatterSolver < Solver
    methods (Access = protected)
        function [vpos, vdat, xform, time] = get_solver_input(obj)
            [vpos, ~, ~, time] = get_solver_input@Solver(obj);             % call superclass method to get positions and time of adcp data
          
            % get backscatter data
            vdat = obj.adcp.backscatter;

            % get transformation matrix
            xform = ones([1, size(vdat,2), size(vdat,3)]);
        end
    end
end