classdef InstrumentMatrixProvider < InstrumentMatrixProvider
    methods(Access = protected)
        function val = get_beam_orientation_matrix(obj, adcp)
            val = obj.i2b_matrix(adcp);
            val(:,:,:,4) = [];
            val = -val;
            % Note that the instrument matrix is defined for the ADCP in a
            % downlooking position, which means that the positive z
            % component results in a vector pointing toward the adcp.
            % Negating it, we have a vector point away from the ADCP.
        end 
        function val = get_vertical_range_to_cell(~,adcp)
            % rdi defines the instrument matrix for a downward looking ADCP
            % so the vertial range is increasing downwards
            val = -(adcp.distmidfirstcell+reshape(0:max(adcp.ncells)-1,[],1).*adcp.cellsize);
        end
    end
end