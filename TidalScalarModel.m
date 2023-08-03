classdef TidalScalarModel < TidalModel & ScalarModel

    methods(Access = public)
        function M = get_model(obj, d_time)
            M = get_model@TidalModel(obj, d_time);
        end
    end
    methods(Access = protected)
        function val = get_npars(obj)
            val = get_npars@TidalModel(obj);
        end

        function names = get_names(obj)
            names = get_names@TidalModel(obj);
        end
    end
end