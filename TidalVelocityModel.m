classdef TidalVelocityModel < TidalModel & VelocityModel
% Tidal velocity model
    methods
        function varargout = get_model(varargin)
            varargout = cell(nargout,1);
            [varargout{:}] = get_model@TidalModel(varargin{:});
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
