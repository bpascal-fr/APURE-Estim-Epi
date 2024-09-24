% Evaluate the linear memory functions of a nonstatinary autoregressive
% process on a vector of observations over the entire observation window:
%   
%    Psi_Y(t) = sum_s Psi(s) * Y(t-s) for s = 1, ..., tau,
%
% for Equation (8) of Pascal & Vaiter (2024).
%
% Border effects at t = 1 are handled by cropping and normalizing the
% coefficients of the memory functions as described in:
%
% - Du, J., Pascal, B., & Abry, P. (2023, August). Compared performance of 
% Covid19 reproduction number estimators based on realistic synthetic data. 
% GRETSI’23 XXIXème Colloque Francophone De Traitement Du Signal Et Des Images.
%
% Default coefficients of the memory functions chosen as the daily 
% discretized serial interval function of COVID-19 estimated in:
%
% - D. Cereda, M. Tirani, F. Rovida, V. Demicheli, M. Ajelli, P. Poletti, 
% F. Trentini, G. Guzzetta, V. Marziano, A. Barone et al. (2020). The early
% phase of the COVID-19 outbreak in Lombardy, Italy. 
% Preprint arXiv:2003.09320
%
% - F. Riccardo, M. Ajelli, X. D. Andrianou, A. Bella, M. Del Manso, 
% M. Fabiani, S. Bellino, S. Boros, A. M. Urdiales, V. Marziano et al. (2020).
% Epidemiological characteristics of COVID-19 cases and estimates of the 
% reproductive numbers 1 month into the epidemic, Italy, 28 January 
% to 31 March 2020. Euro Surveillance.
%
% General nonstationary autoregressive model introduced in:
%
% - Pascal, B., Vaiter, S. (2024, September). Risk Estimate under a
% Nonstationary Autoregressive Model for Data-Driven Reproduction Number
% Estimation. Preprint. arXiv:2409.14937.

% Implementation B. Pascal, J. Du and P. Abry
% March 2024


function [Psi_Y, Y, Psi] = Psi_normal(Y,Psi)

    % Inputs:  - Y: observations of size 1 x T+1
    %          - Psi: coefficients of the linear memory functions (optional) 
    %            (default Gamma distribution of mean 6.6 days and standard
    %            deviation 3.5 days cropped at 25 days and normalized
    %            modeling the daily discretized serial interval function of COVID-19)
    %
    % Outputs: - Y_Psi: memory term (not defined on first time step) hence of size 1 x T
    %          - Y: observations with first time step cropped, hence of size 1 x T
    %          - Psi: coefficients of the linear memory functions
    %          
    

    % Handle Y possibly of size T x 1
    [d1,d2] = size(Y);
    if d2 == 1
        Y = reshape(Y,1,d1);
    end

    if nargin <=1

        % Coefficients of the linear memory function mimicking the daily discretized
        % serial interval function of COVID-19 modeled by a Gamma distribution with
        % - mean: 6.6 days
        % - standard deviation 3.5 days
        % truncated at 25 days and normalized.

        shape    = 1/0.28;
        scale    = 1.87;
        tau      = 25;
        Psi      = gampdf(0:tau,shape,scale);
        Psi      = Psi/sum(Psi); % normalize the weights applied to past tau_phi infection counts

    else

        tau  = length(Psi)-1; % -1 because one zero has been added to Phi for correct conv.

    end


    % compute the convolution of Y by Psi
    Psi_Y             = convn(Y,Psi);
    Psi_Y             = Psi_Y(:,1:size(Y,2));

    % recompute the first tau_phi elements with normalized Phi
    Psi              = reshape(Psi,length(Psi),1);
    Psi_Y(1:tau)  = 0;
    for t = 2:tau

        Phi_crop    = Psi(1:t)./sum(Psi(1:t));
        fY          = fliplr(Y(:,1:t));
        Psi_Y(:,t)  = fY*Phi_crop;

    end

    % Crop the first day at which Psi_Y = 0
    Y               = Y(:,2:end);
    Psi_Y           = Psi_Y(:,2:end);

    % Return Y, Psi_Y of size T x 1 if input Y is univariate of size T+1 x 1
    if d2 == 1
        Y     = reshape(Y,d1-1,d2);
        Psi_Y = reshape(Psi_Y,d1-1,d2);
    end
end