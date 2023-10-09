classdef ExternalDataSolver < Solver
    properties
        position(:,3) double  % x,y,z position of the data
        time(:,1) datetime  % time the data were measured
        data(:,1) double  % data values
        xform(:,:) double  % transformation matrix to transform the data
        water_level_object(1,1) % water levels
    end
    methods(Access = protected)
        function [pos, dat, xform, time, wl] = get_solver_input(obj)
            pos = obj.position;
            dat = obj.data;
            xform = obj.xform;
            time = obj.time;
            wl = obj.water_level_object.get_water_level(obj.time);
        end
    end
end