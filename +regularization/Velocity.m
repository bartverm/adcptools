classdef Velocity <...
        regularization.Regularization &...
        matlab.mixin.Heterogeneous
% Generic class for velocity regularization.
%
%   All classes derived for regularization.Velocity can be combined to
%   weakly enforce several regularizations.
%   
%   methods:
%   get_all_regs - return an object of all velocity regression classes
%
%   methods (protected):
%   mustBeVelocityModel - throws an error if model is not a velocity model
%
%   Available velocity models:
%   regularization.InternalContinuity - enforce continuity within cell
%   regularization.ExternalContinuity - enforce continuity between cells
%   regularization.VelocityCoherence - enforce equal values between cells
%   regularization.VelocityConsistency - enforce equal derivatives
%   regularization.Kinematic - enforce kinematic boundary conditions
%
%   class inherits properties and methods from
%   regularization.Regularization
%
%   see also: regularization.Regularization, 
%     regularization.InternalContinuity, regularization.ExternalContinuity,
%     regularization.VelocityCoherence, regularization.VelocityConsistency,
%     regularization.Kinematic 
    methods(Static)
        function regs = get_all_regs(varargin)
        % returns all available velocity regularizations
        %
        %   regs = regularization.Velocity.get_all_regs(...) returns alls
        %   velocity regularizations and passes given arguments to the
        %   velocity regularization classes for construction.
        %
        %   see also: regularization.Velocity,
        %       regularization.Regularization
            regs = {regularization.InternalContinuity(varargin{:}),...
                    regularization.ExternalContinuity(varargin{:}),...
                    regularization.VelocityCoherence(varargin{:}),...
                    regularization.VelocityConsistency(varargin{:}),...
                    regularization.Kinematic(varargin{:})};
            regs = regularization.Regularization.reorganize_regs(regs);
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