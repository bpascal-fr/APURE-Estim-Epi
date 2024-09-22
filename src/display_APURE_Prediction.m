% Display the oracles associated to the prediction risk and the estimated
% reproduction coefficients obtained using the regularization parameters
% minimizing these oracles.
%
% B. Pascal
% September, 2024

function display_APURE_Prediction(X, lambda, oracle, FontSize)

    % Inputs:  - X: estimated piecewise linear reproduction coefficients
    %            - GT: ground truh (if provided in M.X)
    %            - MLE: estimated maximum likelihood reproduction coefficient
    %            - RISK: reaching lowest true prediction risk (if ground truth is provided in M.X)
    %            - APURE: reaching lowest APURE unbiased prediction risk estimate
    %            - Dates: abstract dates in datetime format or time indices for display (optional, by default 1 to T)
    %            - Name: Prediction or Estimation
    %          - lambda: optimal regularization parameters
    %            - APURE: minimizing APURE unbiased prediction risk estimate
    %            - RISK: minimizing true prediction risk (if ground truth is provided initially)
    %            - GRID: explored range of regularization parameters
    %          - oracles: structure containing for each explored lambda
    %            - APURE: Robustified Finite Difference Monte Carlo APURE unbiased prediction risk estimate
    %            - RISK: Quadratic prediction risk (if ground truth is provided initially)
    %

    if nargin < 4
        FontSize = 22.5 ;
    end

    if ~isfield(X,'Dates');    X.Dates = 1:length(X.APURE);    end

    % Customized colors
    royal  = [0.2549019607843137, 0.4117647058823529, 0.8823529411764706] ;  
    cloud  = [0.65,0.65,0.65] ;
    brick  = [0.6980392156862745, 0.13333333333333333, 0.13333333333333333];
    darko  = [1.0, 0.5490196078431373, 0.0] ;

    

    % Display the oracles
    f31  = figure(31); clf
    if isfield(oracle,'RISK')
        p1 = plot(lambda.GRID,oracle.RISK,'linewidth',2,'color',brick) ;
        hold on
        p2 = plot(lambda.RISK,oracle.RISK(lambda.GRID == lambda.RISK),'.','markersize',35,'Color',brick) ;
        xline(lambda.RISK,'--','LineWidth',2,'Color',brick)
    end
    grid on
    hold on
    fill([lambda.GRID fliplr(lambda.GRID)],[oracle.APURE-oracle.CR fliplr(oracle.APURE-oracle.CR)],darko,'EdgeColor',darko,'LineStyle','--','FaceAlpha',0.3) ;
    p3 = plot(lambda.GRID,oracle.APURE,'linewidth',2,'color',darko) ;
    p4 = plot(lambda.APURE,oracle.APURE(lambda.GRID == lambda.APURE),'.','markersize',35,'Color',darko) ;
    xline(lambda.APURE,'--','LineWidth',2,'Color',darko)
    if isfield(oracle,'RISK')
        l1 = legend([p1,p2,p3,p4],'$\mathcal{P}^\circ$','$\lambda_{\mathcal{P}^\circ}$','$\overline{\mathrm{APURE}_{\mathbf{\zeta}}^{\mathcal{P}}}^N$','$\lambda_{\overline{\mathrm{APURE}_{\mathbf{\zeta}}^{\mathcal{P}}}^N}$','location','best'); 
        l1.Interpreter = 'Latex';
        l1.Position = [0.7592 0.2031 0.1460 0.4384] ;
    else
        l1 = legend([p3,p4],'$\overline{\mathrm{APURE}_{\mathbf{\zeta}}^{\mathcal{P}}}^N$','$\lambda_{\overline{\mathrm{APURE}_{\mathbf{\zeta}}^{\mathcal{P}}}^N}$','location','best'); 
        l1.Interpreter = 'Latex';
        l1.Position    = [0.7593 0.2040 0.1460 0.2632] ;
    end
    l1.FontSize        = FontSize;
    xlabel('$\lambda$','interpreter','latex')
    title('Oracles for the prediction risk and associated optimal $\lambda$','Interpreter','Latex')
    set(gca,'FontSize',FontSize,'ticklabelinterpreter','Latex','XScale','log')
    f31.Position       = [141 329 1033 314];

    % Display the estimates and ground truth if available
    f32  = figure(32); clf

    % Reference vector full of ones
    plot(X.Dates,ones(size(X.APURE)),'-k')
    hold on
    grid on
    Q    = [] ;
    iEst = 1 ; 
    M    = 1 ;

    % Display ground truth if provided
    if isfield(X,'GT')
        q            = plot(X.Dates,X.GT,'linewidth',2,'color',royal) ;
        Q            = [Q, q] ;
        L{iEst}      = '$\overline{\mathrm{X}}_t$' ;
        M            = max(M, max(X.GT)) ;
        iEst         = iEst + 1 ;
    end

    % Display Maximum Likelihood Estimate
    if isfield(X,'MLE')
        q       = plot(X.Dates,X.MLE,'linewidth',1,'color',cloud) ;
        Q       = [Q, q] ;
        L{iEst} = '$\mathrm{X}_t^{\mathrm{ML}}$' ;
        iEst    = iEst + 1 ;
    end

    % Display the estimate obtained minimizing the true risk if available
    if isfield(oracle,'RISK')
        q            = plot(X.Dates,X.RISK,'linewidth',2,'color',brick) ;
        Q            = [Q, q] ;
        L{iEst}      = '$\widehat{\mathrm{X}}_t(\mathrm{Y};\lambda_{\mathcal{P}^\circ})$' ;
        M            = max(M, max(X.RISK)) ;
        iEst         = iEst + 1 ;
    end

    % Display the estimate obtained minimizing the APURE oracle
    q                = plot(X.Dates,X.APURE,'linewidth',2,'color',darko) ;
    Q                = [Q, q] ;
    L{iEst}          = '$\widehat{\mathrm{X}}_t(\mathrm{Y};\lambda_{\overline{\mathrm{APURE}_{\mathbf{\zeta}}^{\mathcal{P}}}^N})$' ;
    M                = max(M, max(X.APURE)) ;
    iEst             = iEst + 1 ;

    % Limits of the y-axis
    ylim([0, 1.05 * M])

    % Display the legend and titles
    l2               = legend(Q,L) ;
    l2.Interpreter   = 'Latex' ;
    if length(Q) == 4
        l2.Position  = [0.7101 0.4589 0.1948 0.4263] ;
    elseif length(Q) == 2
        l2.Position  =  [0.7103 0.6470 0.1948 0.2384] ;
    end
    l2.FontSize      = FontSize ;
    title('Estimates of the reproduction coefficients for optimal parameters $\lambda$','Interpreter','Latex')
    set(gca,'FontSize',FontSize,'ticklabelinterpreter','Latex')
    f32.Position     = [141 329 1033 314];

end