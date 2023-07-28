classdef TidalVelocityModel < TidalModel & VelocityModel
% Tidal velocity model
    methods
        function M = get_model(varargin)
            Mv = get_model@VelocityModel(varargin{:});
            Mt = get_model@TidalModel(varargin{:});
            M = Mv.*Mt;
        end
    end
    methods(Access = protected)
        function val = get_names(obj)
            val = get_names@TidalModel(obj);
        end
        function val = get_npars(obj)
            val = get_npars@TidalModel(obj);
        end
        function val = get_ncomponents(obj)
            val = get_ncomponents@VelocityModel(obj);
        end
    end
end
