% Estimation of the reproduction coefficient under the nonstationary
% autoregressive Poisson model proposed in Equation (11) of Pascal & Vaiter (2024)
%
%         Y(t) | Y(1), ..., Y(t-1) = Poisson(X(t) * Psi_Y(t))
%
% with Psi_Y(t) = sum_s psi(s) * Y(t-s) s = 1, ..., tau, with the sequence
% (psi(s), s = 1, tau) containing the coefficients of the linear memory
% functions.
%
% Maximum Likelihood Estimator writes
%
%         X(t) = Y(t)/Psi_Y(t),
%
% with by convention X(t) = 0 if Psi_Y(t) = 0.
%
% References:
%
% - Pascal, B., Vaiter, S. (2024, September). Risk Estimate under a
% Nonstationary Autoregressive Model for Data-Driven Reproduction Number
% Estimation. Preprint. arXiv:2409.14937.
%
% B. Pascal and S. Vaiter, September 2024.


function X_MLE = X_MaxLikelihood(Y,Psi_Y)

    % Inputs:  - Y: observations under the non
    %          - Psi_Y: memory functions evaluated in the time series realization
    %
    % Output:  - X: estimated maximum likelihood reproduction coefficient

    X_MLE = Y./Psi_Y;

    % Handle trivial estimates
    X_MLE(Psi_Y == 0) = 0;


end