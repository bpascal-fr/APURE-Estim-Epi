% Derivative w.r.t. observations Y of the proximity operator of the 
% Kullback-Leibler divergence between Y and X * Psi_Y at X defined as
%
%      argmin_Q 1/2 * ||X - Q||^2 + gamma * DKL(Y | Q * Psi_Y)
%
% whose closed form expression reads
%
%      prox_gamma*DKL (X) = (X - gamma.*Psi_Y + sqrt((X-gamma.*Psi_Y).^2 + 4*gamma.*Y))/2.
%
% Implementation N. Pustelnik, B. Pascal, C.-G. Lucas and P. Abry
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

function dyprox = dyprox_KLw(X,Y,dY,Psi_Y,dPsi_Y,gamma)

    % Inputs:  - X: reproduction coefficient of size 1 x T
    %          - Y: observations under the nonstationary autoregressive Poisson model of Pascal & Vaiter (2024)
    %          - dY: vector on which the differential is applied
    %          - Psi_Y: linear memory functions evaluted in observations Y
    %          - dPsi_Y: differential of the memory function w.r.t. Y evaluated on dY
    %          - gamma: parameter of the proximity operator (descent step)  (gamma > 0)
    %
    % Outputs: - dyprox: differential of the proximal operator of parameter gamma of 
    %                    the Kullback-Leibler divergence between Y and X * Psi_Y
    %                    evaluated on the vector dY

    num          = sqrt(abs(X-gamma.*Psi_Y).^2 + 4*gamma.*Y);
    dyprox_y     = gamma.*dY./num;
    dyprox_Psi_Y = -gamma/2.*dPsi_Y - gamma/2.*(X-gamma.*Psi_Y).*dPsi_Y./num;
    dyprox       = (dyprox_y + dyprox_Psi_Y);
    
    % Forces output to zero if Y(t) = Psi_Y(t) = 0
    index = find((Psi_Y==0)&(Y==0)) ;

    if ~isempty(index) 

        dyprox(index) = zeros(size(index));

    end

end



