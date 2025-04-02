% Estimation of the reproduction coefficient of a nonstationary
% autoregressive Poisson model with linear memory functions using a
% variational approach consisting in the minimization of the nonsmooth 
% convex functional of Equation (20) of Pascal & Vaiter (2025), in which
% the penalization has been designed to enforce piecewise linearity.
%
% The regularization parameter is chosen manually.
%
% References:
%
% - Abry, P., Pustelnik, N., Roux, S., Jensen, P., Flandrin, P.,
% Gribonval, R., Lucas, C.-G., Guichard, É., Borgnat, P.,
% & Garnier, N. (2020). Spatial and temporal regularization to estimate
% COVID-19 reproduction number R(t): Promoting piecewise smoothness via
% convex optimization. PlosOne, 15(8), e0237901
%
% - Pascal, B., Vaiter, S. (2025). Risk Estimate under a
% Nonstationary Autoregressive Model for Data-Driven Reproduction Number
% Estimation. Preprint. arXiv:2409.14937.
%
% B. Pascal and S. Vaiter, April 2025.


function [X,obj,incr,op] = X_Penalized(Y,Psi_Y,lambda)


    % Minimization of the Poisson penalized log-likelood
    %
    %       DKL(Y | X Psi_Y) + lambda * || D2 X ||_1
    %
    % where DKL stands for the Kullback-Leibler divergence, D2 is the discrete
    % Laplacian operator, || . ||_1 the ell_1-norm defined as the sum of
    % absolute values and lambda > 0 is a regularization parameter.
    %
    % The data fidelity term accounts for the nonstationary autoregressive
    % Poisson model of Equation (11) of Pascal & Vaiter (2025), while the 
    % regularization term enforces a smooth, piecewise-linear behavior 
    % of the reproduction coefficient.
    %
    % Inputs:  - Y: nonstationary autoregressive Poisson observations
    %          - Psi_Y: memory functions evaluated in the observations Y
    %          - lambda: manually tuned regularization parameter
    %
    %
    % Outputs: - X: estimated piecewise linear reproduction coefficients
    %          - obj: values of the objective function w.r.t iterations
    %          - incr: normalized (smoothed) increments w.r.t iterations
    %          - op: linear direct and adjoint operators involved in the regularization term
    

    
    %% RESIZE INPUT 

    [d1,d2]     = size(Y);

    if min(d1,d2) == 1
        
        Y       = reshape(Y,1,max(d1,d2));
        Psi_Y    = reshape(Psi_Y,1,max(d1,d2));

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

    % Regularization term
    objective.regularization = @(y,tau) tau*sum(abs(y(:)));
    prox.regularization      = @(y,tau) prox_L1(y,tau);

    % Linear operator
    filter_def               = 'laplacian' ;
    computation              = 'direct';
    param.type               = '1D' ;
    param.op                 = 'laplacian' ;
    param.lambda             = lambda ;
    op.direct                = @(x)opL(x, filter_def, computation, param) ;
    op.adjoint               = @(x)opLadj(x, filter_def, computation, param) ;
    params.normL             = param.lambda^2 ;

    % Minimization algorithm
    params.dxi               = zeros(size(Z)) ;
    params.xi                = ones(size(Z)) ;
    params.iter              = 1e6 ;
    params.incr              = 'var' ;
    params.prec              = 1e-7 ;
    params.stop              = 'LimSup' ;
    params.win               = 500 ;
    

    %% RUN THE ALGORITHM AND PREPARE OUTPUTS

    % initialization of the primal-dual algorithm
    params.xi         = ones(size(Z)) ;

    % Minimization of the functional with Chambolle-Pock algorithm
    [X, obj, incr]    = PD_ChambollePock_Poisson(Z, objective, op, prox, params);

    % Linear operator involved in the regularization
    param.lambda      = 1;
    op.direct         = @(x)opL(x, filter_def, computation, param);
    op.adjoint        = @(x)opLadj(x, filter_def, computation, param);

    % Handle trivial estimates
    for c = 1:size(Y,1)
        if sum(isnan(Y(c,:))) == size(Y,2)
            X(c,:)    = 0;
        end
    end

    % Resize the output to fit input size
    X                 = reshape(X,d1,d2);

end