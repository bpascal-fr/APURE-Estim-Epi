% Kullback-Leibler divergence between Y and X * Psi_Y
%
%      DKL(Y | X Psi_Y) = sum dkl(Y(t) | X(t) Psi_Y(t))
%
% with
%   dkl(y | p) = y log(y/p) + p - y  if y>0 & p>0,
%   dkl(y | p) = p                   if y = 0 & p >=0
%   dkl(y | p) = Inf                 otherwise.
%
% Implementation N. Pustelnik, CNRS, ENS Lyon
% April 2020
%
% Updated and augmented by P. Abry and B. Pascal
% March 2024
%
% Adpated to general nonstationary autoregressive model by B. Pascal
% September 2024
%
% Reference:
%
% - Pascal, B., Vaiter, S. (2024, September). Risk Estimate under a
% Nonstationary Autoregressive Model for Data-Driven Reproduction Number
% Estimation. Preprint. arXiv:.


function KLw = KLw(X,Y,Psi_Y)

    % Inputs:  - X: reproduction coefficient of size 1 x T
    %          - Y: observations under the nonstationary autoregressive Poisson model of Pascal & Vaiter (2024)
    %          - Psi_Y: linear memory functions evaluted in observations Y
    %
    % Outputs: - DKL: Kullback-Leibler divergence between Y and X * Psi_Y

    % Find indices for which dkl is finite
    j = find(Y>0  & Psi_Y.*X>0 );
    k = find(Y==0 & X>=0);

    % Compute the sum of individual dkl 
    if length(j) + length(k) < length(X)

        KLw  = Inf;

    else

        DKLj      = Psi_Y(j) .* X(j) - Y(j) + Y(j).*log(Y(j)./(Psi_Y(j).*X(j)));
        DKLk      = Psi_Y(k) .* X(k);

        KLw  = sum(DKLj(:))+sum(DKLk(:));

    end

end