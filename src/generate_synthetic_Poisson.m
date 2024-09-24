% Generate synthetic nonstationary autoregressive Poisson time series with
% linear memory functions following the Model (11) of Pascal & Vaiter (2024)
%
%      Y(t) ~ alpha * Poisson(X(t) * Psi(t) / alpha)
%
% where Psi(t) = sum psi(s) * Y_(t-s), s = 1, ..., tau, with the sequence
% (psi(s), s = 1, tau) containing the coefficients of the linear memory
% functions. In this example, psi mimics the discretized serial interval
% distribution of COVID-19, modeled as a Gamma distribution of
% mean 6.6 days and standard deviation 3.5 days cropped at tau = 25 days.
%
%
% References:
%
% - Du, J., Pascal, B., & Abry, P. (2023). Performances comparées
% d’estimateurs du coefficient de reproduction de la Covid19 à l’aide de
% données synthétiques réalistes. GRETSI’23 XXIXème Colloque Francophone De
% Traitement Du Signal Et Des Images.
%
% - Pascal, B., Vaiter, S. (2024, September). Risk Estimate under a
% Nonstationary Autoregressive Model for Data-Driven Reproduction Number
% Estimation. Preprint. arXiv:2409.14937.
%
% B. Pascal and S. Vaiter, September 2024.


function [Y, Psi_Y, M] = generate_synthetic_Poisson(X, sett)

    % Inputs: - X: ground truth reproduction coefficient vector of size 1 x T (trivial: ones(1,T))
    %         - sett: parameters of the model, structure containing
    %                   sett.Y0: initial number (default: 12)
    %                   sett.alpha: scale parameter of the Poisson model (default: 1, corresponding to usual Poisson distribution)
    %                   sett.FontSize: font size in the plots (default FontSize = 22.5)
    %                   sett.Dates: abstract dates in datetime format or time indices for display (optional, by default 1 to T)
    %                   sett.Psi: coefficients of the linear memory functions (default: daily discretized Covid19 serial interval function)
    %
    %
    % Outputs: - Y: synthetic nonstationary autoregressive Poisson time series with linear memory functions
    %          - Psi_Y: memory functions evaluated in the time series realization
    %          - M: model parameters, structure containing
    %                   M.X: ground truth reproduction coefficient
    %                   M.Y0: initial observation
    %                   M.Y_SY: synthetic observations including the initial observation
    %                   M.Psi: coefficients of the linear memory functions
    %                   M.alpha: scale parameter of the Poisson model
    %                   M.Dates: abstract dates in datetime format for display

    if nargin < 2

        % intial observation
        Y0      = 12;

        % scale parameter
        alpha    = 1;

        % plot font size
        FontSize = 22.5;

        % temporal axis
        Dates     = 1:length(X);

    else

        if ~isfield(sett,'Y0');       sett.Y0 = 12;                end
        if ~isfield(sett,'alpha');    sett.alpha = 1;              end
        if ~isfield(sett,'FontSize'); sett.FontSize = 22.5;        end
        if ~isfield(sett,'Dates');    sett.Dates = 1:length(X);    end

        % intial observation
        Y0      = sett.Y0;

        % scale parameter
        alpha    = sett.alpha;

        % plot font size
        FontSize = sett.FontSize;

        % temporal axis
        Dates    = sett.Dates;

    end

    % Coefficients of the linear memory functions mimicking
    if ~isfield(sett,'Psi') % default: Covid19 serial interval function
        tau     = 25;                        % memory horizon
        shape   = 1/0.28;                    % shape parameter
        scale   = 1.87;                      % scale parameter
        Psi     = gampdf(0:tau,shape,scale); % discretized Gamma probability density function
    else
        if ~(sett.Psi(1) == 0)
            tau = length(sett.Psi);
            Psi = reshape(sett.Psi,1,tau);
            Psi = [0, Psi];
        else
            tau = length(sett.Psi) - 1;
            Psi = reshape(sett.Psi,1,tau+1);
        end
    end
    Psi     = Psi/sum(Psi);             % normalize the memory function (optional)
    fPsi    = fliplr(Psi(2:end));

    % Prepare synthetic time series
    T       = length(X);
    Y_SY    = zeros(1,T+1);
    Y_SY(1) = Y0;
    Psi_Y   = zeros(1,T);

    % Threshold to enforce positive Poisson intensity
    Thr     = 3;
    
    % Generate synthetic observations recursively
    for t = 2:T+1

        if t < tau+2
            tPsi      = Psi(2:t)./sum(Psi(2:t));
            tY        = fliplr(Y_SY(1:t-1));
            Psi_Y(t-1) = sum(tY.*tPsi);
        else
            Psi_Y(t-1) = sum(Y_SY(t-tau:t-1).*fPsi);
        end

        P_sy    = X(t-1)*Psi_Y(t-1);
        P_sy    = max(P_sy,Thr)/alpha;
        Y_SY(t) = alpha*random('Poisson',P_sy);

    end

    % Remove initial observation
    Y           = Y_SY(2:end);

    % Store model parameters
    M.X        = X ;
    M.Y0       = Y0 ;
    M.Y_SY     = Y_SY ;
    M.Psi      = Psi ;
    M.alpha    = alpha ;
    M.Dates    = Dates ;

    %% DISPLAY SYNTHETIC COUNTS

    f2              = figure(2); clf
    p               = plot(Dates, Y,'-','linewidth',2,'color','black') ;
    grid on ; hold on
    q               = plot(Dates, Psi_Y,'-.','linewidth',2,'color','black') ;
    leg             = legend([p,q],'$\mathrm{Y}_t$','$\Psi_t(Y)$','location','best');
    leg.Interpreter = 'Latex';
    leg.Color       = 'none';
    leg.FontSize    = FontSize;
    VX              = [Dates(1) Dates(end)] ;
    xlim(VX)
    title(strcat("Nonstationary autoregressive Poisson time series, $\alpha = ",num2str(alpha,5),"$"),'Interpreter','Latex')
    ylabel('observations','Interpreter','Latex')
    set(gca,'ticklabelinterpreter','Latex','fontsize',FontSize,'color','None')
    f2.Position     = [141 329 1033 314];

end