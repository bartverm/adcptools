classdef SideLobeFilter < Filter
% Filter to exclude velocity affected by sidelobes and below the bed
    methods
        function obj=SideLobeFilter()
            obj.description='Side lobe filter';
        end
    end
    methods(Access=protected)
        function bad=bad_int(obj,adcp)
            validateattributes(adcp,{'VMADCP'},{'scalar'})
            if isscalar(obj)
               bed_range=-permute(max(adcp.bed_offset,[],3),[1 2 4 3]);            % get highest bed level detected at each ensemble
               bed_range=bed_range(:,:,3)-adcp.cellsize/2;                          % remove half cell size such that the entire cell is not affected by side lobes
               tm=-adcp.xform(CoordinateSystem.Beam, CoordinateSystem.Earth);   % get beam orientation
               tm(:,:,:,[1 2 4])=[];                                                % keep vertical component of beam orientation
               min_lev=bed_range.*tm;                                                 % compute minimul elevation velocity needs to be valid
               vel_pos=adcp.depth_cell_offset(CoordinateSystem.Earth);            % get velocity position
               vel_pos=vel_pos(:,:,:,3);
               bad=vel_pos<min_lev;                                                 % mark bad cells
            else
                bad=bad_int@Filter(obj,adcp);
            end
        end
    end
end