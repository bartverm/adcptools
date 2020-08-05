classdef TimeBasedVelocitySolver < VelocitySolver
    properties

    end
    methods
        function obj=TimeBasedVelocitySolver(varargin)
            obj=obj@VelocitySolver(varargin{:});
            for cnt_arg=1:nargin
                cur_arg=varargin{cnt_arg};
                if isa(cur_arg, 'VMADCP') || isa(cur_arg, 'Mesh') || isa(cur_arg,'Bathymetry') || isa(cur_arg, 'Filter') || isa(cur_arg, 'XSection')
                    continue
                else
                    warning(['Unhandled input of type: ', class(cur_arg), ' on construction of TimeBasedVelocitySolver object'])
                end
            end
        end
        function vel=get_velocity(obj)
            vpos=obj.adcp.depth_cell_position;
            filters=[obj.filter; obj.adcp.filters];
            vpos(filters.bad(obj.adcp))=nan;
            vpos=nanmean(vpos,3);
            [~, n_pos]=obj.xs.xy2sn(vpos(:,:,:,1),vpos(:,:,:,2));
            zb_pos=obj.bathy.get_bed_elev(vpos(:,:,:,1), vpos(:,:,:,2));
            sig_pos=1-vpos(:,:,:,3)./zb_pos;
            idx=repmat(obj.mesh.index(n_pos, sig_pos), 1, 1, 3);
            fcomp=cumsum(ones(size(idx)),3);
            vel=nan(obj.mesh.ncells, 3);
            vel_data = obj.adcp.water_velocity(CoordinateSystem.Earth);
            for cc = 1:obj.mesh.ncells
                for comp=1:3
                    vel(cc,comp) = nanmean(vel_data(idx==cc & fcomp==comp)); 
                end
            end
        end
    end
end