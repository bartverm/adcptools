classdef SolverOptions < handle
    %SOLVER_OPTIONS Contains input parameters for solving for mesh-based
    % parameters

    properties
        % Solver options

        algorithm = 'pcg'; %lscov, pcg
        % Algorithm that performs the actual inversion of the linear
        % system(s) of equations. 'lscov' reverts back to the solver of
        % Vermeulen et al. 2014. 'pcg' (preconditioned conjugate gradient
        % method) assembles a global system of equations and is suitable
        % for regularized inversion.

        % in case lscov
        lscov_opts = [];

        % in case pcg
        preconditioner_opts (1,1) struct = struct('michol', 'on', 'type', 'ict', 'droptol', 1e-3,'diagcomp', 1e-3);
        set_diagcomp (1,1) double = 1;
        % Dynamically sets the modified incomplete Cholesky preconditioner
        % according to: alpha = max(sum(abs(A),2)./diag(A))-2;
        %               preconditioner_opts.diagcomp = max(sum(abs(A),2)./diag(A))-2;
        % This is done at the point of calling pcg

        pcg_tol (1,1) double {mustBePositive} = 1e-9;
        pcg_iter (1,1) double {mustBeInteger} = 1000;

        % Initialization of iterative solver using base vector p, which is
        % a solution to the regularized problem with regularization
        % parameters reg_pars(1)
        use_p0 (1,1) double = 1;

        % Cross-validation options

        cv_mode = 'random'; %none, random, omit_cells, omit_time
        % Generalized cross-validation mode. 'none': no cross-validation.
        % 'random', 'omit_cells', 'omit_time': Perform data splitting in
        % training and validation sets, according to the exclusion rules
        % provided by the three parameters hereunder:

        % in case of cv_mode = 'random'
        training_perc = .90; % percentage of data in training set
        cv_iter = 1; % Re-partitioning iterations into training and validation sets (affects also gen_analysis)

        % in case of cv_mode = 'omit_cells' or 'omit_time'
        omit_cells (1,:) double {mustBeInteger} = []; % in case of omit_cells: indices of the cells to omit
        omit_time (1,:) double = []; % in case of omit_time: time interval to omit in the training set


        % Regularization options

        reg_pars (1,:) cell = {[100, 100, 5, 5, 100]}; % Regularization parameters. Can be cell array of double vectors of same length
        
        % Affects both reg_pars and reg_pars_sens
        no_flow = {'right', 'surface', 'left', 'bottom'}; % Declare no-flow boundaries (affects all computations)

        % Sensitivity study: Generalization error
        gen_analysis = 0;
        reg_vary = 'coupled'; % 'coupled', 'full'
        % If set to 'coupled': Variation of lambda_c and lambda_d. If set
        % to 'full': independent variation of all regularization hyperparameters (not
        % recommended)
        reg_relative_weights = [1, 1, 1, 1, 1];
        force_zero (1,5) double = [0,0,0,0,0]; % Set index i to one if you want to set lambda_i to zero

        reg_iter (1,5) double = [4, 4, 4, 4, 4]; % iterations per regularization parameter.
        % Only reg_iter(1) and reg_iter(3) will play a role if reg_vary =
        % 'coupled'.
        res_near_zero (1,5) double {mustBePositive} = .05*ones(1,5); % The smaller this number,
        % the more the regularization par will be clustered
        % towards zero in the sensitivity experiment

        max_reg_pars (1,5) double {mustBePositive} = [1e3, 1e3, 1e3, 1e3, 1e3]; % vector of maximal regularization hyperparameters
        % in the sensitivity experiment
        min_reg_pars (1,5) double = [0, 0, 0, 0, 0]; % Fixed

        reg_pars_sens (1,5) cell = {nan, nan, nan, nan, nan};

        % TO INCLUDE: SENSITIVITY TO ADDITIVE NOISE options.
    end

    methods
        function obj = SolverOptions(varargin)
            %SolverOptions

            % Overwrite default solver options
            for ia = 1:2:nargin
                obj.(varargin{ia}) = varargin{ia+1};
            end

            % prepare generalization error analysis
            for i = 1:5
                ymax = helpers.symlog(obj.max_reg_pars(i), obj.res_near_zero(i));
                ymin = helpers.symlog(obj.min_reg_pars(i), obj.res_near_zero(i));
                obj.reg_pars_sens{i} = helpers.symexp(linspace(ymin, ymax, obj.reg_iter(i)), obj.res_near_zero(i))';
            end
            
        end

        function reg_pars_sens_vec = vectorize_reg_pars(obj)
            if strcmp(obj.reg_vary, 'coupled')
                L = obj.reg_pars_sens(1, [1,3]); % Only vary two reg. parameters
                n = length(L);
                [L{:}] = ndgrid(L{end:-1:1});
                L = cat(n+1,L{:});
                L = fliplr(reshape(L,[],n));
                RP = [L(:, 1), L(:,1), L(:,2), L(:,2), L(:,1)]; % Coupling of parameters
            elseif strcmp(obj.reg_vary, 'full')
                L = obj.reg_pars_sens; % Vary all reg. parameters
                n = length(L);
                [L{:}] = ndgrid(L{end:-1:1});
                L = cat(n+1,L{:});
                RP = fliplr(reshape(L,[],n));
            else
                error('Invalid reg_vary option')
            end

            for i = 1:size(RP,2)
                RP(:,i) = obj.reg_relative_weights(i).*RP(:,i);
                if obj.force_zero(i)    
                    RP(:,i) = 0;
                end
            end


            reg_pars_sens_vec = RP;
        end
    end
end

