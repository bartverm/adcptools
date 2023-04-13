classdef TaylorBasedRegularization < Regularization
    methods(Access = protected)
        function assemble_matrix_private(obj,varargin)
            if ~obj.model_is_taylor
                warning("TaylorBasedRegularizon:NoTaylorModel",...
                    "Regularization requires a TaylorModel");
            end
        end
        function tf = model_is_taylor(obj)
            tf = isa(obj.model,"TaylorModel");
        end
    end
end