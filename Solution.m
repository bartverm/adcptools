classdef Solution < handle & helpers.ArraySupport
    %MODELPARAMETERS Class for capturing model parameters as output of
    %get_parameters

    %   Detailed explanation goes here

    properties(SetAccess = ?Solver)
        M (:,:) double % matrix M such that b = Mp

        b (:,1) double % rhs of system of eqs (b = Mp)

        p (:,:) double % model parameters (b = Mp)

        mesh (1,1) Mesh = SigmaZetaMesh

        regularization (1,:) regularization.Regularization

        model (1,1) DataModel = TaylorScalarModel

        opts (1,1) SolverOptions

        ns (:,:) double

    end
    properties(Dependent)
        pars
    end
    methods
        function obj = Solution(varargin)
            obj = obj@helpers.ArraySupport(varargin{:})
        end
        function pars = get.pars(obj)
            [np, nsols] = size(obj.p);
            ncells = obj.mesh.ncells;
            npars = np/ncells;
            pars = reshape(obj.p, npars, ncells, nsols);
            pars = permute(pars, [2 1 3]);
        end
        function varargout = get_data(obj)
            varargout = cell(1,nargout);
            [varargout{:}] = obj.model.get_data(obj.pars);
        end
    end
end
