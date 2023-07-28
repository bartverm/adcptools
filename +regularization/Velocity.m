classdef Velocity <...
        regularization.Regularization &...
        matlab.mixin.Heterogeneous
    methods
        function obj = Velocity(varargin)
        end
    end
    methods(Static)
        function regs = get_all_regs(varargin)
            regs = {regularization.InternalContinuity(varargin{:}),...
                    regularization.ExternalContinuity(varargin{:}),...
                    regularization.VelocityCoherence(varargin{:}),...
                    regularization.VelocityConsistency(varargin{:}),...
                    regularization.Kinematic(varargin{:})};
            siz = cellfun(@size, regs, 'UniformOutput',false);
            assert(isequal(siz{:}),...
                'Size of generated regularizations should be equal')
            siz = siz{1};
            if isequal(siz, [1 1])
                regs = [regs{:}];
                return
            end
            out = cell(siz);
            for co = 1:numel(out)
                tmp = cellfun(@(x) x(co), regs, 'UniformOutput',false);
                tmp = [tmp{:}];
                out{co} = tmp;
            end
            regs = out;
        end
    end
    methods(Access = protected)
        function mustBeVelocityModel(obj)
            val = all(cellfun(@(x) isa(x,'VelocityModel'), {obj.model}),...
                'all');
            assert(val, "Regularization requires a VelocityModel");
        end
    end
end