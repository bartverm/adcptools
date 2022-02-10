classdef InstrumentMatrixProvider < InstrumentMatrixProvider
    methods(Static, Access = protected)
        function val = get_cart2beam_to_orient(tm)
            val = tm;
            val(:,:,:,4) = [];
            val = -val;
        end 
    end
end