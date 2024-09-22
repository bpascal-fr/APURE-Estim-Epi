% Estimation of the reproduction coefficient of a nonstationary
% autoregressive Poisson model with linear memory functions with a
% piecewise linearity prior through minimization of the nonsmooth convex
% regularizing functional of Equation (20) of Pascal & Vaiter (2024).
%
% The regularization parameter is selected automatically by the
% minimization of the (robustfied) Finite Difference Monte Carlo APURE
% unbiased prediction risk estimate of Equations (28) and (36) of
% Pascal & Vaiter (2024).
%
%
%
% References:
%
% - Abry, P., Pustelnik, N., Roux, S., Jensen, P., Flandrin, P.,
% Gribonval, R., Lucas, C.-G., Guichard, É., Borgnat, P.,
% & Garnier, N. (2020). Spatial and temporal regularization to estimate
% COVID-19 reproduction number R(t): Promoting piecewise smoothness via
% convex optimization. PlosOne, 15(8), e0237901
%
% - Pascal, B., Vaiter, S. (2024, September). Risk Estimate under a
% Nonstationary Autoregressive Model for Data-Driven Reproduction Number
% Estimation. Preprint. arXiv:.
%
% B. Pascal and S. Vaiter, September 2024.


function [X, lambda, oracle] = APURE_Prediction(Y,Psi_Y,M,opts)


    % Minimization of the Poisson penalized log-likelood
    %
    %       DKL(Y | X Psi_Y) + lambda * || D2 X ||_1
    %
    % where DKL stands for the Kullback-Leibler divergence, D2 is the discrete
    % Laplacian operator, || . ||_1 the ell_1-norm defined as the sum of
    % absolute values and lambda > 0 is a regularization parameter.
    %
    % The data fidelity term accounts for the nonstationary autoregressive
    % Poisson model, while the regularization term enforces a smooth,
    % piecewise-linear behavior of the reproduction coefficient.
    %
    % The parameter lambda is selected by the minimization of the
    % robustfied Finite Difference Monte Carlo APURE unbiased prediction
    % risk estimate w.r.t. lambda:
    %
    %      lambda_APURE = argmin_lambda 1/N sum_n APURE_P[zeta_n](Y; lambda | alpha)
    %
    % where zeta_1, ..., zeta_N are i.i.d. standard Gaussian random vectors.
    % The minimization is performed on a grid of L values of lambda from
    % lambda_min to lambda_max.
    %
    % Inputs:  - Y: nonstationary autoregressive Poisson observations
    %          - Psi_Y: memory functions evaluated in the observations Y
    %          - M: structure containing parameters of the model
    %            - X: ground truth reproduction coefficient (optional, if not provided the true error is not computed)
    %            - Y0: initial observation (optional, if not provided Y(1) is used instead)
    %            - Psi: coefficients of the linear memory functions (optional, if not provided mimic COVID-19 serial interval function)
    %            - alpha: scale parameter of the Poisson model (optional, if not provided use 1 corresponding to standard Poisson)
    %            - Dates: abstract dates in datetime format or time indices for display (optional, by default 1 to T)
    %          - opts: structure containing the model knowledge, the
    %          properties of the unbiased risk estimate and of the grid
    %          search on lambda
    %            - N: number of Monte Carlo vectors used in the robustified APURE estimate (default: 10)
    %            - L: number of values of lambda explored (default: 60)
    %            - lambda_min: smallest value of lambda explored (default: 1e-2)
    %            - lambda_min: largest value of lambda explored (default: 1e4)
    %            - flag: if 'none' no progression bar (optional)
    %
    %
    % Outputs: - X: estimated piecewise linear reproduction coefficients
    %            - GT: ground truh (if provided in M.X)
    %            - MLE: estimated maximum likelihood reproduction coefficient
    %            - RISK: reaching lowest true prediction risk (if ground truth is provided in M.X)
    %            - APURE: reaching lowest APURE unbiased prediction risk estimate
    %            - Dates: abstract dates in datetime format for display (if provided in M.Dates)
    %            - Name: Prediction or Estimation
    %          - lambda: optimal regularization parameters
    %            - APURE: minimizing APURE unbiased prediction risk estimate
    %            - RISK: minimizing true prediction risk (if ground truth is provided in M.X)
    %            - GRID: explored range of regularization parameters
    %          - oracles: structure containing for each explored lambda
    %            - APURE: Robustified Finite Difference Monte Carlo APURE unbiased prediction risk estimate
    %            - RISK: Quadratic prediction risk (if ground truth is provided in M.X)


    %% DEFAULTS OPTIONS

    if nargin == 3
        opts     = struct;
    end

    % Parameters of the robustified Finite Difference Monte Carlo unbiased risk estimate
    if ~isfield(opts,'N'),          opts.N          = 10 ;   end

    % Parameters of the grid search on the regularization parameter lambda
    if ~isfield(opts,'L'),          opts.L          = 60 ;   end
    if ~isfield(opts,'lambda_min'), opts.lambda_min = 1e-2 ; end
    if ~isfield(opts,'lambda_max'), opts.lambda_max = 1e4 ;  end

    % Name of the oracle for displaying waiting bar
    if isfield(opts,'flag')
        if strcmp(opts.flag,'none')
            opts = rmfield(opts,'flag');
        end
    else
        opts.flag = 'APURE in Prediction';
    end

    %% RESIZE INPUT

    [d1,d2]     = size(Y);

    if min(d1,d2) == 1

        Y       = reshape(Y,1,max(d1,d2));
        Psi_Y   = reshape(Psi_Y,1,max(d1,d2));

    end

    %% NORMALIZE OBSERVATIONS AND MEMORY TERMS

    scale       = std(Y,[],2);   % scale of infection counts
    Z           = Y./scale;
    Psi_Z       = Psi_Y./scale;

    %% PARAMETERS OF THE FUNCTIONAL AND OF THE PROXIMAL ALGORITHM

    params      = struct ;

    % Regularizing functional

    % Data-fidelity term
    params.mu                = 0;
    objective.fidelity       = @(y,Y) KLw(y,Y,Psi_Z);
    prox.fidelity            = @(y,Y,tau) prox_KLw(y,Y,Psi_Z,tau);
    dprox.fidelity           = @(x,dx,y,gamma) dprox_KLw(x,dx,y,Psi_Z,gamma);

    % Regularization term
    objective.regularization = @(y,tau) tau*sum(abs(y(:)));
    prox.regularization      = @(y,tau) prox_L1(y,tau);
    dprox.regularization     = @(y,dy,tau) dprox_L1(y,dy,tau);

    % Linear operator
    filter_def               = 'laplacian' ;
    computation              = 'direct';
    param.type               = '1D' ;
    param.op                 = 'laplacian' ;

    % Minimization algorithm
    params.dxi               = zeros(size(Z)) ;
    params.xi                = ones(size(Z)) ;
    params.iter              = 1e6 ;
    params.incr              = 'var' ;
    params.prec              = 1e-7 ;
    params.stop              = 'LimSup' ;
    params.win               = 500 ;

    %% GRID SEARCH TO FIND THE MINIMUM OF THE APURE UNBIASED PREDICTION RISK ESTIMATE

    % Data-driven estimate of the Poisson scale parameter
    if ~isfield(M,'alpha'), M.alpha = 0.1 * std(Y) ; end

    % Design the regularization paremeters grid to be explored
    Lambdas      = logspace(log10(opts.lambda_min),log10(opts.lambda_max),opts.L) ;

    % Constant term of APURE (independent of lambda and of Monte Carlo vectors)
    C            = sum(Y.^2) - M.alpha*sum(Y) ;

    % Initialize oracles
    APURE_MC     = zeros(opts.N,opts.L) ;
    if isfield(M,'X')
        RISK     = zeros(1,opts.L) ;
    end

    % Initialize progression bar
    if isfield(opts,'flag')
        f        = waitbar(0,['Monte Carlo vector: 0 / ',num2str(opts.N),'; Parameter \lambda: 0 / ',num2str(opts.L)],'Name',['Minimization of ',opts.flag]);
    end

    % Explore the grid for N realizations of the Monte Carlo vector
    for n = 1 : opts.N

        % Sample the standard Gaussian Monte Carlo vector
        dY_SY                = randn(size(M.Y_SY)) ;
        if isfield(M,'Psi')
            [dPsi_Y,dY]      = Psi_normal(dY_SY, M.Psi) ;
            dPsi_Z           = dPsi_Y / scale ;
            dZ               = dY / scale ;
        else
            [dPsi_Y,dY]      = Psi_normal(dY_SY) ;
            dPsi_Z           = dPsi_Y / scale ;
            dZ               = dY / scale ;
        end

        % Define the differentiated proximity operator depending on dPsi_Y
        dprox.yfidelity  = @(x,y,dy,gamma) dyprox_KLw(x,y,dy,Psi_Z,dPsi_Z,gamma);

        for ell = 1 : opts.L

            if sum(isnan(Z))      == length(Z)

                % Handle trivial estimates
                X_lambda          = zeros(size(Z)) ;
                dX_lambda         = zeros(size(Z)) ;

            else

                % Minimization of the functional with Chambolle-Pock algorithm
                param.lambda          = Lambdas(ell) ;
                op.direct             = @(x)opL(x, filter_def, computation, param) ;
                op.adjoint            = @(x)opLadj(x, filter_def, computation, param) ;
                params.normL          = param.lambda^2 ;
                [X_lambda, dX_lambda] = dPD_ChambollePock_Poisson(Z, dZ, objective, op, prox, dprox, params);

            end

            % Compute the Finite Difference Monte Carlo APURE unbiased prediction risk estimate
            APURE_MC(n,ell)       =  norm(Psi_Y .* X_lambda , 'fro')^2 - 2 * sum(Y .* Psi_Y .* X_lambda) + 2 * M.alpha * sum(Y .* Psi_Y .* dX_lambda .* dY) + C ;

            % Compute true prediction risk if ground truth available
            if n == 1
                if isfield(M,'X')
                    RISK(ell)     = norm(Psi_Y .* M.X - Psi_Y .* X_lambda,'fro')^2 ;
                end
            end

            % Show progression of the algorithm
            if isfield(opts,'flag')
                waitbar(((n-1) * opts.L + ell) / (opts.N * opts.L),f,['Monte Carlo vector: ',num2str(n),' / ',num2str(opts.N),'; Parameter \lambda: ',num2str(ell),' / ',num2str(opts.L)]);
            end

        end

    end

    % End of the exploration of the grid, close waiting bar
    pause(1)
    close(f)
    
    % Robustified oracle
    APURE        = mean(APURE_MC,1) ;

    % Gaussian confidence region across Monte Carlo vectors
    CR           = 1.96 / sqrt(opts.N) * std(APURE_MC,[],1) ;

    % Regularization parameter minimizing robustfied APURE
    opt_ell      = find(APURE == min(APURE)) ;
    lbd_opt      = Lambdas(opt_ell(1)) ;

    if isfield(M,'X')
        % Best regularization parameter in terms of prediction risk
        best_ell = find(RISK == min(RISK)) ;
        lbd      = Lambdas(best_ell(1)) ;
    end

    % Optimal estimates
    if sum(isnan(Z))      == length(Z)
        X_opt             = zeros(size(Z));
    else
        disp('Computing the estimate with optimal parameter minimizing APURE ...')
        param.lambda      = lbd_opt ;
        op.direct         = @(x)opL(x, filter_def, computation, param) ;
        op.adjoint        = @(x)opLadj(x, filter_def, computation, param) ;
        params.normL      = param.lambda^2 ;
        X_opt             = PD_ChambollePock_Poisson(Z, objective, op, prox, params);
        clc
    end
    if isfield(M,'X')
        if sum(isnan(Z))      == length(Z)
            X_best            = zeros(size(Z));
        else
            disp('Computing the estimate with optimal parameter minimizing the true risk ...')
            param.lambda      = lbd ;
            op.direct         = @(x)opL(x, filter_def, computation, param) ;
            op.adjoint        = @(x)opLadj(x, filter_def, computation, param) ;
            params.normL      = param.lambda^2 ;
            X_best            = PD_ChambollePock_Poisson(Z, objective, op, prox, params);
            clc
        end
    end

    % Store estimates and ground truth if available
    if isfield(M,'X')
        X.GT              = M.X ;
    end
    X.MLE                 = X_MaxLikelihood(Y, Psi_Y) ; 
    if isfield(M,'X')
        X.RISK            = X_best ;
    end
    X.APURE               = X_opt ;
    if isfield(M,'Dates')
        X.Dates           = M.Dates ;
    end
    X.Name                = 'Prediction' ;

    % Store optimal regularization parameters and whole explored grid
     lambda.APURE         = lbd_opt ;
     if isfield(M,'X')
        lambda.RISK       = lbd ;
    end
    lambda.GRID           = Lambdas ;

    % Store oracles values and confidence region at all point of the grid
    oracle.APURE    = APURE ;
    oracle.CR       = CR ;
    if isfield(M,'X')
        oracle.RISK = RISK ;
    end


end