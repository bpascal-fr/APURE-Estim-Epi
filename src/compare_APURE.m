% Display the oracles associated to the estimation risk and the estimated
% reproduction coefficients obtained using the regularization parameters
% minimizing these oracles.
%
% B. Pascal
% September, 2024

function compare_APURE(X1,X2,FontSize)

    % Inputs:  - X1, X2: estimated piecewise linear reproduction coefficients
    %            - APURE: reaching lowest APURE unbiased estimation risk estimate
    %            - RISK: reaching lowest true estimation risk (if ground truth is provided in initially)
    %            - GT: ground truh (if provided initially)
    %            - Dates: abstract dates in datetime format or time indices for display (optional, by default 1 to T)
    %            - Name: Prediction or Estimation

    if nargin < 3
        FontSize = 22.5 ;
    end

    if ~isfield(X1,'Dates');    X1.Dates = 1:length(X1.APURE);    end
    if ~isfield(X2,'Dates');    X2.Dates = 1:length(X2.APURE);    end

    % Customized colors
    royal    = [0.2549019607843137, 0.4117647058823529, 0.8823529411764706] ;
    cloud    = [0.65,0.65,0.65] ;
    brick    = [0.6980392156862745, 0.13333333333333333, 0.13333333333333333];
    darko    = [1.0, 0.5490196078431373, 0.0] ;
    olive    = [0.4196078431372549, 0.5568627450980392, 0.13725490196078433];
    orchid   = [0.8549019607843137, 0.4392156862745098, 0.8392156862745098];


    if ~ (length(X1.APURE) == length(X2.APURE))
        error('Not same time period: comparison makes no sense.')
    else 
        if isfield(X1,'GT') && isfield(X2, 'GT')
            if ~ (sum(X1.GT == X2.GT) == length(X1.GT))
                error('Not same ground truth: comparison makes no sense.')
            end
        end
    end

    

    % Display the estimates and ground truth if available
    f5       = figure(5); clf

    % Reference vector full of ones
    plot(X1.Dates,ones(size(X1.APURE)),'-k')
    hold on
    grid on
    Q        = [];
    iEst     = 1 ;
    M        = 1 ;

    % Display ground truth if provided
    if isfield(X1,'GT')
        q            = plot(X1.Dates,X1.GT,'linewidth',2,'color',royal) ;
        Q            = [Q, q] ;
        L{iEst}      = '$\overline{\mathrm{X}}_t$' ;
        iEst         = iEst + 1 ;
        M            = max(M, max(X1.GT)) ;
    else
        if isfield(X2,'GT')
            q        = plot(X2.Dates,X2.GT,'linewidth',2,'color',royal) ;
            Q        = [Q, q] ;
            L{iEst}  = '$\overline{\mathrm{X}}_t$' ;
            iEst     = iEst + 1 ;
            M            = max(M, max(X2.GT)) ;
        end
    end

    % Display Maximum Likelihood Estimate
    if isfield(X1,'MLE')
        q            = plot(X1.Dates,X1.MLE,'linewidth',2,'color',cloud) ;
        Q            = [Q, q] ;
        L{iEst}      = '$\mathrm{X}_t^{\mathrm{ML}}$' ;
        iEst         = iEst + 1 ;
    else
        if isfield(X2,'GT')
            q        = plot(X2.Dates,X2.MLE,'linewidth',2,'color',cloud) ;
            Q        = [Q, q] ;
            L{iEst}  = '$\mathrm{X}_t^{\mathrm{ML}}$' ;
            iEst     = iEst + 1 ;
        end
    end

    % Display the first estimates
    if strcmp(X1.Name,'Prediction')
        if isfield(X1,'RISK')
            q        = plot(X1.Dates,X1.RISK,'linewidth',2,'color',brick) ;
            Q        = [Q, q] ;
            L{iEst}  = '$\widehat{\mathrm{X}}(\mathrm{Y};\lambda_{\mathcal{P}^\circ})$' ;
            iEst     = iEst + 1 ;
            M        = max(M, max(X1.RISK)) ;
        end
        q            = plot(X1.Dates,X1.APURE,'linewidth',2,'color',darko) ;
        Q            = [Q, q] ;
        L{iEst}      = '$\widehat{\mathrm{X}}_t(\mathrm{Y};\lambda_{\overline{\mathrm{APURE}_{\mathbf{\zeta}}^{\mathcal{P}}}^N})$' ;
        iEst         = iEst + 1 ;
        M            = max(M, max(X1.APURE)) ;
    elseif strcmp(X1.Name,'Estimation')
        if isfield(X1,'RISK')
            q        = plot(X1.Dates,X1.RISK,'linewidth',2,'color',olive) ;
            Q        = [Q, q] ;
            L{iEst}  = '$\widehat{\mathrm{X}}_t(\mathrm{Y};\lambda_{\mathcal{E}^\circ})$' ;
            iEst     = iEst + 1 ;
            M        = max(M, max(X1.RISK)) ;
        end
        q            = plot(X1.Dates,X1.APURE,'linewidth',2,'color',orchid) ;
        Q            = [Q, q] ;
        L{iEst}      = '$\widehat{\mathrm{X}}_t(\mathrm{Y};\lambda_{\overline{\mathrm{APURE}_{\mathbf{\zeta}}^{\mathcal{E}}}^N})$' ;
        iEst         = iEst + 1 ;
        M            = max(M, max(X1.APURE)) ;
    else
        warning('The first estimates have not been identified and will hence be ignored.')
    end

    % Display the second estimates
    if strcmp(X2.Name,'Prediction')
        if isfield(X2,'RISK')
            q        = plot(X2.Dates,X2.RISK,'linewidth',2,'color',brick) ;
            Q        = [Q, q] ;
            L{iEst}  = '$\widehat{\mathrm{X}}(\mathrm{Y};\lambda_{\mathcal{P}^\circ})$' ;
            iEst     = iEst + 1 ;
            M        = max(M, max(X2.RISK)) ;
        end
        q            = plot(X2.Dates,X2.APURE,'linewidth',2,'color',darko) ;
        Q            = [Q, q] ;
        L{iEst}      = '$\widehat{\mathrm{X}}_t(\mathrm{Y};\lambda_{\overline{\mathrm{APURE}_{\mathbf{\zeta}}^{\mathcal{P}}}^N})$' ;
        iEst         = iEst + 1 ;
        M            = max(M, max(X2.APURE)) ;
    elseif strcmp(X2.Name,'Estimation')
        if isfield(X2,'RISK')
            q        = plot(X2.Dates,X2.RISK,'linewidth',2,'color',olive) ;
            Q        = [Q, q] ;
            L{iEst}  = '$\widehat{\mathrm{X}}_t(\mathrm{Y};\lambda_{\mathcal{E}^\circ})$' ;
            iEst     = iEst + 1 ;
            M        = max(M, max(X2.RISK)) ;
        end
        q            = plot(X2.Dates,X2.APURE,'linewidth',2,'color',orchid) ;
        Q            = [Q, q] ;
        L{iEst}      = '$\widehat{\mathrm{X}}_t(\mathrm{Y};\lambda_{\overline{\mathrm{APURE}_{\mathbf{\zeta}}^{\mathcal{E}}}^N})$' ;
        iEst         = iEst + 1 ;
        M            = max(M, max(X2.APURE)) ;
    else
        warning('The second estimates have not been identified and will hence be ignored.')
    end

    % Limits of the y-axis
    ylim([0, 1.1 * M])
    
    % Display the legend and titles
    leg              = legend(Q,L) ;
    leg.Interpreter  = 'Latex' ;
    leg.FontSize     = FontSize ;
    title('Estimates of the reproduction coefficients for optimal parameters $\lambda$','Interpreter','Latex')
    set(gca,'FontSize',FontSize,'ticklabelinterpreter','Latex')
    f5.Position      = [141 329 1033 314];
    if length(Q)     == 6
        leg.Position = [0.7103 0.2222 0.1948 0.6631] ;
    elseif length(Q) == 3
        leg.Position = [0.7103 0.5075 0.1948 0.3779] ;
    end

end