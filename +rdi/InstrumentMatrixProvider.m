classdef InstrumentMatrixProvider < InstrumentMatrixProvider
    methods(Access = protected)
        function val = get_beam_orientation_matrix(obj, adcp)
            val = obj.i2b_matrix(adcp);
            val(:,:,:,4) = [];
            val = -val;
        end 
    end
end