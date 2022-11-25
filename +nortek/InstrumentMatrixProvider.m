classdef InstrumentMatrixProvider < InstrumentMatrixProvider
% Nortek defines the instrument to beam matrix for an upward looking ADCP.
% Beam velocity is positive away from ADCP.
    methods(Access = protected)
        function val = get_beam_orientation_matrix(obj, adcp,varargin)
            val = obj.get_i2b_matrix(adcp);
            val = cat(4, val(:,:,:,1:2),...
                    cat(3,...
                    val(:,:,1,3), ...
                    val(:,:,2,4), ...
                    val(:,:,3,3), ...
                    val(:,:,4,4)));
        end 
        function val = get_vertical_range_to_cell(~,adcp)
            val = (adcp.distmidfirstcell+...
                reshape(0:max(adcp.ncells)-1,[],1).*adcp.cellsize);
        end
    end
end