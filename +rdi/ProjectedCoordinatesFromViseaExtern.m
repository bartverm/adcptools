classdef ProjectedCoordinatesFromViseaExtern < ProjectedCoordinatesProvider
    properties
        description='';
    end
    methods(Access=protected)
        function tf=get_has_data(~,adcp)
            tf=isfield(adcp.raw,'VISEA_Extern') && all(isfield(adcp.raw.VISEA_Extern,{'Northing','Easting'}));
        end
        function pos=get_horizontal_position(~,adcp)
            x=adcp.raw.VISEA_Extern.Easting;
            y=adcp.raw.VISEA_Extern.Northing;
            pos=[x;y];
        end
    end
end