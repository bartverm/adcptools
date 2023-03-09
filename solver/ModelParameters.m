classdef ModelParameters < handle
    %MODELPARAMETERS Class for capturing model parameters as output of
    %get_parameters
    
    %   Detailed explanation goes here
    
    properties
        M (:,:) double % matrix M such that b = Mp

        b (:,1) double; %rhs of system of eqs (b = Mp)

        p (:,:) double; % model parameters (b = Mp)

        reg (1,1) Regularization

        opts (1,1) SolverOptions
    end
    
    methods
        function obj = ModelParameters(varargin)
            % Overwrite default options
            for ia = 1:2:nargin
                obj.(varargin{ia}) = varargin{ia+1};
            end
        end
       
        % HERE, INCLUDE VARIOUS POST-PROCESSING FUNCTIONS.

    end
end

