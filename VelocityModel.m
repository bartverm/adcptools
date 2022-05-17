classdef VelocityModel < handle
% Model for Cartesian velocity components
%
%   This base class implements a dummy model, i.e. each velocity is
%   represented by the mean. Subclasses can implement linear models for the
%   velocity. They should implement a linear mapping between model
%   parameters and cartesian velocity components. Subclasses can
%   reimplement the get_model methods and the get_npars methods to
%   implement a new velocity model.
%
%   VelocityModel properties:
%   npars - number of parameters in the velocity model
%
%   VelocityModel methods:
%   get_velocity - return velocity based on model parameters
%
%  see also: TaylorExpandedVelocity, VelocitySolver

    properties(Dependent, SetAccess=private)
        % VelocityModel/npars (read only) number of parameteris in the model
        %
        %   3x1 row vector holding the number of model parameters for each
        %   velocity components. For VelocityModel this always returns [1 1
        %   1] since this class implements the simplest model, i.e. the
        %   velocity is equal to the mean velocity.
        %
        %   see also: VelocityModel, get_model
        npars
    end
    methods

        function [vel, cov_vel] = get_velocity(obj, pars, cov_pars, d_time, d_s, d_n, d_z, d_sigma)
        % Compute velocity from model parameters.
        %
        %   vel = get_velocity(obj, pars) computes velocity based on model
        %   parameters. 
        %
        %   [vel, cov_vel] = get_velocity(obj, pars, cov_pars) also compute
        %   covariance matrix of velocity based on covariance of model
        %   parameters.
        %   
        %   [vel, cov_vel] = get_velocity(obj, pars, cov_pars, d_time, d_s, 
        %       d_n, d_z, d_sigma) optionally specify the coordinate
        %       offsets at which velocities will be computed. These default
        %       to zeros if not given. They are all column vectors with
        %       same number of rows as the rows in pars.
        %
        %   see also: VelocitySolver   

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
            [Mu, Mv, Mw] = obj.get_model(d_time, d_s, d_n, d_z, d_sigma);
            
            % reconstruct model matrix with kron product
            np = obj.npars;
            nin = size(pars,1);
            z1 = zeros(nin, np(1));
            z2 = zeros(nin, np(2));
            z3 = zeros(nin, np(3));
            M = cat(3,...
                cat(2, Mu, z2, z3),...
                cat(2, z1, Mv, z3),...
                cat(2, z1, z2, Mw));

            % apply model to obtain vel (permutations for use of pagemtimes)
            M = permute(M,[3, 2 ,1]);
            pars = permute(pars, [2, 3, 1]);
            vel = pagemtimes(M,pars);

            % apply model to obtain covariance matrix
            cov_pars = permute(cov_pars, [2, 3, 1]);
            cov_vel = pagemtimes(cov_pars,'none', M,'transpose');
            cov_vel = pagemtimes(M, cov_vel);

            % transform output back to format of inputs
            vel = ipermute(vel, [2, 3, 1]);
            cov_vel = ipermute(cov_vel, [2, 3, 1]);
        end
        
        function val=get.npars(obj)
            val=obj.get_npars();
        end

        function [Mu, Mv, Mw] = get_model(~, d_time, ~, ~, ~, ~)
        % build velocity model matrices
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
        %   For VelocityModel, the matrices hold a column of ones, meaning
        %   that each measured velocity component is an estimate of the
        %   final mean velocity.
        %   Reimplementing this in a subclass allows to define other, more
        %   comples, but linear velocity models.
        %
        %   d_time, d_s, d_n, d_z, d_sigma are the distances of the
        %   measured velocity data from the center of the cell, or the time
        %   difference with the mesh time.
        %
        %   see also: VelocitySolver, TaylorExpandedVelocity.

            [Mu, Mv, Mw]=deal(ones(numel(d_time),1));
        end
    end
    methods(Access=protected)
        function val=get_npars(~)
        % return number of parameters as a 3x1 row vector, with the number
        % of parameters for the x,y and z components, respectively.
            val=[1 1 1];
        end
    end
end