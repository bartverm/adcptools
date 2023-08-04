classdef TaylorVelocityModel < VelocityModel & TaylorModel
% Taylor expanded velocity model
    methods
        function M = get_model(obj,varargin)
            M = obj.get_model@TaylorModel(varargin{:});
        end
    end
    methods(Access = protected)
        function val = get_names(obj)
            val = get_names@TaylorModel(obj);
        end
        function val = get_npars(obj)
            val = get_npars@TaylorModel(obj);
        end
        function val = get_ncomponents(obj)
            val = get_ncomponents@VelocityModel(obj);
        end
    end
end