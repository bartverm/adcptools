classdef DataModel < handle
% Model for Cartesian velocity components
%
%   This base class implements a dummy model, i.e. each velocity is
%   represented by the mean. Subclasses can implement linear models for the
%   velocity. They should implement a linear mapping between model
%   parameters and cartesian velocity components. Subclasses can
%   reimplement the get_model methods and the get_npars methods to
%   implement a new velocity model.
%
%   DataModel properties:
%   npars - number of parameters in the velocity model
%
%   DataModel methods:
%   get_velocity - return velocity based on model parameters
%
%  see also: TaylorModel, TidaModel, DataSolver

    properties(Dependent, SetAccess=private)
        % DataModel/npars (read only) number of parameteris in the model
        %
        %   1xNCOMP row vector holding the number of model parameters for each
        %   data components. For DataModel this always returns [1 1
        %   1] since this class implements the simplest model, i.e. the
        %   velocity is equal to the mean velocity.
        %
        %   see also: DataModel, get_model
        npars

        % DataModel/ncomponents (read only) components in data
        %
        %   3x1 row vector holding the number of model parameters for each
        %   velocity components. For DataModel this always returns [1 1
        %   1] since this class implements the simplest model, i.e. the
        %   velocity is equal to the mean velocity.
        %
        %   see also: DataModel, get_model
        ncomponents        
    end
    methods

        function [dat, cov_dat] = get_data(obj, pars, cov_pars, d_time, d_s, d_n, d_z, d_sigma)
        % Compute data values from model parameters.
        %
        %   dat = get_data(obj, pars) computes data based on model
        %   parameters. 
        %
        %   [dat, cov_dat] = get_data(obj, pars, cov_pars) also compute
        %   covariance matrix of data based on covariance of model
        %   parameters.
        %   
        %   [dat, cov_dat] = get_velocity(obj, pars, cov_pars, d_time, d_s, 
        %       d_n, d_z, d_sigma) optionally specify the coordinate
        %       offsets at which data will be computed. These default
        %       to zeros if not given. They are all column vectors with
        %       same number of rows as the rows in pars.
        %
        %   see also: DataModel, Solver

            t_var=zeros(size(pars,1),1);
            mult=ones(size(t_var));
            if nargin < 4
                d_time = t_var; 
            else 
                if isscalar(d_time)
                    d_time = repmat(d_time, size(pars, 1), 1);
                end
            end
            if nargin < 5
                d_s = t_var; 
            else 
                d_s = d_s .* mult; 
            end
            if nargin < 6
                d_n = t_var; 
            else 
                d_n = d_n .* mult; 
            end
            if nargin < 7
                d_z = t_var;  
            else 
                d_z = d_z .* mult; 
            end
            if nargin < 8
                d_sigma = t_var;  
            else 
                d_sigma = d_sigma .* mult; 
            end
            M = obj.get_model(d_time, d_s, d_n, d_z, d_sigma);
            
            % reconstruct model matrix with kron product
            np = obj.npars;
            nin = size(pars,1);
            nc = obj.ncomponents;
            Mnew = zeros(nin,nc,sum(np));
            cum_pars = cumsum([0 np]);
            for c_comp = 1:obj.ncomponents
                Mnew(:,c_comp,cum_pars(c_comp)+1:cum_pars(c_comp+1))=...
                M(:,1:np(c_comp),c_comp);
            end
            M=Mnew;

            % apply model to obtain vel (permutations for use of pagemtimes)
            dat = helpers.matmult(M,pars);

            % apply model to obtain covariance matrix
            cov_dat = helpers.matmult(cov_pars, permute(M,[1,3,2]));
            cov_dat = helpers.matmult(M, cov_dat);
        end
        
        function val=get.npars(obj)
            val=obj.get_npars();
        end
        function val=get.ncomponents(obj)
            val=obj.get_ncomponents();
        end

        function M = get_model(~, d_time, ~, ~, ~, ~)
        % build velocity model matrices
        %
        %   M is NEnsembles x NPars x Ncomp
        %
        %   [Mu, Mv, Mw] = get_model(obj, d_time, d_s, d_n, d_z,
        %   d_sigma) given the offsets of the velocity data to the mesh
        %   cell coordinates return the model matrices which describe the
        %   relation between model parameters and cartesian velocity.
        %   For each cartesian velocity component a model is given in Mu,
        %   Mv and Mw for the x, y, and z component respectively. Mu, Mv,
        %   and Mw are matrices with each row corresponding to a measured
        %   velocity component and each column being the known coefficient
        %   for a model parameter.
        %   For DataModel, the matrices hold a column of ones, meaning
        %   that each measured velocity component is an estimate of the
        %   final mean velocity.
        %   Reimplementing this in a subclass allows to define other, more
        %   complex, but linear velocity models.
        %
        %   d_time, d_s, d_n, d_z, d_sigma are the distances of the
        %   measured velocity data from the center of the cell, or the time
        %   difference with the mesh time.
        %
        %   see also: Solver, TaylorModel, TidalModel.

            M = ones(numel(d_time),1,3);
        end
    end
    methods(Access=protected)
        function val=get_npars(~)
        % return number of parameters as a 3x1 row vector, with the number
        % of parameters for the x,y and z components, respectively.
            val=[1 1 1];
        end
        function val = get_ncomponents(~)
            val = 3;
        end
    end
end