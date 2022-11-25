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
                % get highest bed level detected at each ensemble
               bed_range=permute(max(adcp.bed_offset,[],3),[1 2 4 3]);            
                
               % retain only third dimension and negate to get distance
               bed_range=-bed_range(:,:,3);

               % get beam orientation
               beam_or=adcp.beam_orientation_matrix;
               tm = adcp.xform(CoordinateSystem.Earth,...
                   CoordinateSystem.Instrument,...
                   'UseTilts', true); 
               % permute for proper rotation
               beam_or = helpers.matmult(tm(:,:,1:3,1:3),...
                   permute(beam_or,[1,2,4,3])); 
               beam_or = permute(beam_or,[1,2,4,3]);

                % compute minimum acceptable velocity level
               min_lev=bed_range.*beam_or(:,:,:,3)+adcp.cellsize/2;
               vel_pos=adcp.depth_cell_offset(CoordinateSystem.Earth);
               vel_pos=vel_pos(:,:,:,3);
               bad=vel_pos<min_lev;
            else
                bad=bad_int@Filter(obj,adcp);
            end
        end
    end
end