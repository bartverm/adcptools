classdef ModelParameters < handle
    %MODELPARAMETERS Class for capturing model parameters as output of
    %get_parameters
    
    %   Detailed explanation goes here
    
    properties
        model (1,1) DataModel = DataModel; %Underlying model for scalar or vector quantities
        
        M (1,:) double = {}; % matrix M such that b = Mp

        b (:,1) double = []; %rhs of system of eqs (b = Mp)

        p (:,:) double = []; % model parameters (b = Mp)

        C (1,:) cell = {}; % regularization matrices

        cell_idx (:,1) double = [];
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

