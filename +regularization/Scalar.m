classdef Scalar <...
        regularization.Regularization &...
        matlab.mixin.Heterogeneous

    methods(Static)
        function regs = get_all_regs(varargin)
            regs = {regularization.ScalarCoherence(varargin{:}),...
                    regularization.ScalarConsistency(varargin{:})};
            regs = regularization.Regularization.reorganize_regs(regs);
        end
    end
end