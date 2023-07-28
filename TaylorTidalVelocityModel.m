classdef TaylorTidalVelocityModel < TaylorTidalModel & VelocityModel
% Taylor expanded tidal velocity model
    methods
        function M = get_model(varargin)
            Mv = get_model@VelocityModel(varargin{:});
            Mtt = get_model@TaylorTidalModel(varargin{:});
            M = Mv.*Mtt;
        end
    end
    methods(Access = protected)
        function val = get_names(obj)
            val = get_names@TaylorTidalModel(obj);
        end
        function val = get_npars(obj)
            val = get_npars@TaylorTidalModel(obj);
        end
        function val = get_ncomponents(obj)
            val = get_ncomponents@VelocityModel(obj);
        end
    end
end