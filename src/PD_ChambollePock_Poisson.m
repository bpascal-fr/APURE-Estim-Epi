% Primal-dual algorithm from Chambolle and Pock to solve a nonsmooth
% convex optimization problem of the form
%
%    minimize f(x) + r(Lx)
%
% with f,r nonsmooth convex proper semicontinuous functions and L a bounded
% linear operator.
%
% If f is strong-convexity, the accelerated version is implemented.
%
% References:
%
% - Chambolle, A., & Pock, T. (2011). A first-order primal-dual algorithm
% for convex problems with applications to imaging. Journal of mathematical
% imaging and vision, 40, 120-145.
%
% Stopping criterion on increments of R from
% - Du, J., Pascal, B., & Abry, P. (2023, August). Compared performance of
% Covid19 reproduction number estimators based on realistic synthetic data.
% GRETSI’23 XXIXème Colloque Francophone De Traitement Du Signal Et Des Images.
%
% Implementation N. Pustelnik, CNRS, ENS Lyon
% April 2020
%
% Updated and augmented by P. Abry and B. Pascal
% March 2024


function [x, obj, incr] = PD_ChambollePock_Poisson(Y, objective, op, prox, opts)


    % Inputs:  - Y: observed corrupted data
    %          - objective: tools to compute the objective function
    %               fidelity: data fidelity term f
    %               regularization: penalization term g
    %          - op: linear operators
    %               direct: bounded linear operator L
    %               adjoint: adjoint operator of L
    %          - prox: proximity operators
    %               fidelity: proximity operator of the data fidelity term f
    %               regularization: proximity operator of the penalization term g
    %          - opts: structure containing the properties of the regularizing functional and the parameters of the minimization algorithm
    %               normL: squared operator norm of L
    %               mu: strong-convexity modulus of the data fidelity term
    %               xi: initialization of x
    %               iter: maximal number of iterations
    %               incr: 'var' for increments on iterates, 'obj' increments on objective function
    %               prec: tolerance for the stopping criterion
    %               stop: 'LimSup' smoothed increments over win past iterates, or 'Primal' pointwise increments
    %               win: length of smoothing window (500 by default)
    %               flag: name of the estimator to display waiting bars
    %
    %
    % Outputs: - x: minimizer of the objective function
    %          - obj: values of the objective function w.r.t iterations
    %          - incr: normalized smoothed increments w.r.t iterations

    % Primal-dual descent steps
    gamma = 0.99;
    tau   = gamma/sqrt(opts.normL);
    sig   = gamma/sqrt(opts.normL);
    theta = 1;

    % Initializing primal and dual variables
    x  = opts.xi;
    y  = op.direct(x);
    x0 = x;
    bx = x;

    % Criteria of convergence
    obj   = zeros(1,opts.iter);
    tinc  = obj;
    incr  = obj;

    % Initialize the number of iteration and the current precision
    i     = 0;
    incrc = opts.prec + 1;

    % Initialize progression bar
    if isfield(opts,'flag')
        f = waitbar(0,['Iterations: 0 / ',num2str(opts.iter)],'Name',['Computing ',opts.flag]);
    end

    %% Iterations of Chambolle-Pock algorithm
    
    % If all data are zeros then the estimate is zero else optimization can be run
    if sum(Y(:) == 0) == numel(Y)

        x    = zeros(size(opts.xi));
        obj  = obj(1);
        incr = incr(1);

    else

    while (incrc > opts.prec) && (i < opts.iter)

        i = i+1;

        % Store the reproduction number for computing increments
        if strcmp(opts.incr,'var')

            Xm = x(:,1:size(Y,2));

        end

        % Update the primal variable
        tmp = y + sig*op.direct(bx);
        y   = tmp - sig*prox.regularization(tmp/sig, 1/sig);

        % Update the dual variable
        x   = prox.fidelity(x0 - tau * op.adjoint(y),Y,tau);

        % Update the the descent steps if acceleration is possible
        if opts.mu >= 0

            theta = (1+2*opts.mu*tau)^(-1/2);
            tau   = theta*tau;
            sig   = sig/theta;

        end

        % Update dual auxiliary variable
        bx = x + theta*(x - x0);
        x0 = x;

        % Compute the objective function
        obj(i) = objective.fidelity(x,Y) +  objective.regularization(op.direct(x),1);

        % Compute the increments
        if i > 1

            if strcmp(opts.incr,'obj')

                objm          = obj(i-1);
                tinc(i-1)     = abs(obj(i) - objm)/abs(objm);

                if strcmp(opts.stop,'Primal')

                    incr(i-1) = tinc(i-1);

                elseif strcmp(opts.stop,'LimSup')

                    ind_past  = max(1,i-opts.win);
                    incr(i-1) = max(tinc(ind_past:i-1));

                end

            else

                R             = x(:,1:size(Y,2));
                if isempty(find(Xm(:) > 0,1))
                    if i == 2
                        tinc(i-1) = incrc;
                    else
                        tinc(i-1) = tinc(i-2);
                    end
                else
                    tinc(i-1)     = max(abs(Xm(Xm(:) > 0)-R(Xm(:) > 0))./max(1e-2,Xm(Xm(:) > 0)));
                end

                if strcmp(opts.stop,'Primal')

                    incr(i-1) = tinc(i-1);

                elseif strcmp(opts.stop,'LimSup')

                    ind_past  = max(1,i-opts.win);
                    incr(i-1) = max(tinc(ind_past:i-1));

                end

            end

            incrc = incr(i-1);

        end

        % Show progression of the algorithm
        if isfield(opts,'flag')
            if mod(i,1000) == 0
                waitbar(i/opts.iter,f,['Iterations: ',num2str(i),' / ',num2str(opts.iter)]);
            end
        end

    end

    % Show ending of the algorithm
    if isfield(opts,'flag')
        if i == opts.iter
            waitbar(1,f,'Iteration: 100 % (maximal number of iterations reached)');
            pause(1)
            close(f)
            disp([opts.flag,' estimator has reached the maximal number of iterations: ',num2str(opts.iter),' and stoped at precision: ',num2str(incrc,'%1.e'),'.'])
        else
            waitbar(i/opts.iter,f,['Iterations: ',num2str(i),' / ',num2str(opts.iter),' (required precision reached)']);
            pause(1)
            close(f)
            disp([opts.flag,' estimator has reached the required precision: ',num2str(opts.prec,'%1.e'),'.'])
        end
    end

    % Crop objective function and increments at the effective number of iterations
    obj  = obj(1:i);
    incr = incr(1:i-1);

    end

end

