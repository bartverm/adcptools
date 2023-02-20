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
        
        pcg_tol = 1e-9;
        pcg_iter = 1000;
        
        % Initialization of iterative solver using base vector p0
        use_p0 (1,1) double = 1;

        % Cross-validation options

        cv_mode = 'none'; %none, random, omit_cells, omit_time
        % Generalized cross-validation mode. 'none': no cross-validation.
        % 'random', 'omit_cells', 'omit_time': Perform data splitting in
        % training and validation sets, according to the exclusion rules
        % provided by the three parameters hereunder:

        % in case of cv_mode = 'random'
        training_perc = .90; % percentage of data in training set
        cv_iter = 1; % Re-partitioning iterations into training and validation sets

        % in case of cv_mode = 'omit_cells' or 'omit_time'
        omit_cells (1,:) double = []; % in case of omit_cells: indices of the cells to omit
        omit_time (1,2) double = []; % in case of omit_time: time interval to omit in the training set

        % Sensitivity study options

        reg_vary = 'coupled'; % 'coupled', 'full'
        % If set to 'coupled': Variation of lambda_c and lambda_d. If set
        % to 'full': independent variation of all regularization hyperparameters (not
        % recommended)

        reg_iter (1,5) double = [4, nan, 4, nan, nan]; % iterations per regularization parameter. 
        res_near_zero (1,1) double {mustBePositive} = .05; % The smaller this number, 
        % the more the regularization hyperpareters will be clustered
        % towards zero in the hyperparameter sensitivity experiment

        max_reg_pars (1,5) double {mustBePositive} = [1e3, 1e3, 1e3, 1e3, 1e3]; % vector of maximal regularization hyperparameters
        % in the sensitivity experiment
        min_reg_pars (1,5) double = [0, 0, 0, 0, 0]; % Fixed
                
        reg_pars = {100, 100, 5, 5, 100}; % Regularization parameters

        reg_pars_compare = {[0, 100, 1000], [0, 100, 1000], [0, 20, 50], [0, 20, 50], [0, 100, 1000]}; % Cell of
        % vectors of regularization parameters, which will all be used in model inversion.
        reg_pars_sensitivity (1,5) cell = []; % Default reg parameters. Can also be multiple for purposes of comparing.


        force_zero (1,5) double = [0,0,0,0,0]; % Set to one if you want to set lambda_i to zero.

    end
    
    methods
        function obj = SolverOptions(varargin)
            %SOLVER_OPTIONS Construct an instance of this class
            %   Detailed explanation goes here
            %p = inputParser('caseSensitive', 1, )
            ymax = symlog(opts.max_reg_pars, opts.res_near_zero);

            reg_pars = {symexp(linspace(0,ymax,opts.reg_iter(1)), opts.res_near_zero)',...
                symexp(linspace(0,ymax,opts.reg_iter(2)), opts.res_near_zero)',...
                symexp(linspace(0,ymax,opts.reg_iter(3)), opts.res_near_zero)',...
                symexp(linspace(0,ymax,opts.reg_iter(4)), opts.res_near_zero)',...
                symexp(linspace(0,ymax,opts.reg_iter(5)), opts.res_near_zero)'};
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

