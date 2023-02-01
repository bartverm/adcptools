classdef MagneticDeviationInternal < MagneticDeviationModel
    methods
        function val = magnetic_deviation(~,adcp)
            val = -reshape(adcp.raw.Compass.Magnetic_error(:,1),1,[]);
        end
    end
end