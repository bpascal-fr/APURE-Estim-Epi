% Display the oracles associated to the estimation risk and the estimated
% COVID-19 reproduction number obtained using the regularization parameters
% minimizing these oracles.
%
% B. Pascal
% September, 2024

function compare_Covid(R1,R2,FontSize)

    % Inputs:  - R1, R2: estimated piecewise linear COVID-19 reproduction number
    %            - APURE: reaching lowest APURE unbiased estimation risk estimate
    %            - RISK: reaching lowest true estimation risk (if ground truth is provided in initially)
    %            - GT: ground truh (if provided initially)
    %            - Dates: abstract dates in datetime format or time indices for display (optional, by default 1 to T)
    %            - Name: Prediction or Estimation

    if nargin < 3
        FontSize = 22.5 ;
    end

    if ~isfield(R1,'Dates');    R1.Dates = 1:length(R1.APURE);    end
    if ~isfield(R2,'Dates');    R2.Dates = 1:length(R2.APURE);    end

    % Customized colors
    royal    = [0.2549019607843137, 0.4117647058823529, 0.8823529411764706] ;
    cloud    = [0.65,0.65,0.65] ;
    brick    = [0.6980392156862745, 0.13333333333333333, 0.13333333333333333];
    darko    = [1.0, 0.5490196078431373, 0.0] ;
    olive    = [0.4196078431372549, 0.5568627450980392, 0.13725490196078433];
    orchid   = [0.8549019607843137, 0.4392156862745098, 0.8392156862745098];


    if ~ (length(R1.APURE) == length(R2.APURE))
        error('Not same time period: comparison makes no sense.')
    else 
        if isfield(R1,'GT') && isfield(R2, 'GT')
            if ~ (sum(R1.GT == R2.GT) == length(R1.GT))
                error('Not same ground truth: comparison makes no sense.')
            end
        end
    end

    

    % Display the estimates and ground truth if available
    f5       = figure(5); clf

    % Reference vector full of ones
    plot(R1.Dates,ones(size(R1.APURE)),'-k')
    hold on
    grid on
    Q        = [];
    iEst     = 1 ;
    M        = 1 ;

    % Display ground truth if provided
    if isfield(R1,'GT')
        q            = plot(R1.Dates,R1.GT,'linewidth',2,'color',royal) ;
        Q            = [Q, q] ;
        L{iEst}      = '$\overline{\mathrm{R}}_t$' ;
        iEst         = iEst + 1 ;
        M            = max(M, max(R1.GT)) ;
    else
        if isfield(R2,'GT')
            q        = plot(R2.Dates,R2.GT,'linewidth',2,'color',royal) ;
            Q        = [Q, q] ;
            L{iEst}  = '$\overline{\mathrm{R}}_t$' ;
            iEst     = iEst + 1 ;
            M            = max(M, max(R2.GT)) ;
        end
    end

    % Display Maximum Likelihood Estimate
    if isfield(R1,'MLE')
        q            = plot(R1.Dates,R1.MLE,'linewidth',2,'color',cloud) ;
        Q            = [Q, q] ;
        L{iEst}      = '$\mathrm{R}_t^{\mathrm{ML}}$' ;
        iEst         = iEst + 1 ;
    else
        if isfield(R2,'GT')
            q        = plot(R2.Dates,R2.MLE,'linewidth',2,'color',cloud) ;
            Q        = [Q, q] ;
            L{iEst}  = '$\mathrm{R}_t^{\mathrm{ML}}$' ;
            iEst     = iEst + 1 ;
        end
    end

    % Display the first estimates
    if strcmp(R1.Name,'Prediction')
        if isfield(R1,'RISK')
            q        = plot(R1.Dates,R1.RISK,'linewidth',2,'color',brick) ;
            Q        = [Q, q] ;
            L{iEst}  = '$\widehat{\mathrm{R}}(\mathrm{Z};\lambda_{\mathcal{P}^\circ})$' ;
            iEst     = iEst + 1 ;
            M        = max(M, max(R1.RISK)) ;
        end
        q            = plot(R1.Dates,R1.APURE,'linewidth',2,'color',darko) ;
        Q            = [Q, q] ;
        L{iEst}      = '$\widehat{\mathrm{R}}_t(\mathrm{Z};\lambda_{\overline{\mathrm{APURE}_{\mathbf{\zeta}}^{\mathcal{P}}}^N})$' ;
        iEst         = iEst + 1 ;
        M            = max(M, max(R1.APURE)) ;
    elseif strcmp(R1.Name,'Estimation')
        if isfield(R1,'RISK')
            q        = plot(R1.Dates,R1.RISK,'linewidth',2,'color',olive) ;
            Q        = [Q, q] ;
            L{iEst}  = '$\widehat{\mathrm{R}}_t(\mathrm{Z};\lambda_{\mathcal{E}^\circ})$' ;
            iEst     = iEst + 1 ;
            M        = max(M, max(R1.RISK)) ;
        end
        q            = plot(R1.Dates,R1.APURE,'linewidth',2,'color',orchid) ;
        Q            = [Q, q] ;
        L{iEst}      = '$\widehat{\mathrm{R}}_t(\mathrm{Z};\lambda_{\overline{\mathrm{APURE}_{\mathbf{\zeta}}^{\mathcal{E}}}^N})$' ;
        iEst         = iEst + 1 ;
        M            = max(M, max(R1.APURE)) ;
    else
        warning('The first estimates have not been identified and will hence be ignored.')
    end

    % Display the second estimates
    if strcmp(R2.Name,'Prediction')
        if isfield(R2,'RISK')
            q        = plot(R2.Dates,R2.RISK,'linewidth',2,'color',brick) ;
            Q        = [Q, q] ;
            L{iEst}  = '$\widehat{\mathrm{R}}(\mathrm{Z};\lambda_{\mathcal{P}^\circ})$' ;
            iEst     = iEst + 1 ;
            M        = max(M, max(R2.RISK)) ;
        end
        q            = plot(R2.Dates,R2.APURE,'linewidth',2,'color',darko) ;
        Q            = [Q, q] ;
        L{iEst}      = '$\widehat{\mathrm{R}}_t(\mathrm{Z};\lambda_{\overline{\mathrm{APURE}_{\mathbf{\zeta}}^{\mathcal{P}}}^N})$' ;
        iEst         = iEst + 1 ;
        M            = max(M, max(R2.APURE)) ;
    elseif strcmp(R2.Name,'Estimation')
        if isfield(R2,'RISK')
            q        = plot(R2.Dates,R2.RISK,'linewidth',2,'color',olive) ;
            Q        = [Q, q] ;
            L{iEst}  = '$\widehat{\mathrm{R}}_t(\mathrm{Z};\lambda_{\mathcal{E}^\circ})$' ;
            iEst     = iEst + 1 ;
            M        = max(M, max(R2.RISK)) ;
        end
        q            = plot(R2.Dates,R2.APURE,'linewidth',2,'color',orchid) ;
        Q            = [Q, q] ;
        L{iEst}      = '$\widehat{\mathrm{R}}_t(\mathrm{Z};\lambda_{\overline{\mathrm{APURE}_{\mathbf{\zeta}}^{\mathcal{E}}}^N})$' ;
        iEst         = iEst + 1 ;
        M            = max(M, max(R2.APURE)) ;
    else
        warning('The second estimates have not been identified and will hence be ignored.')
    end
    
    % Limits of the x-axis
    xlim([R1.Dates(1), R1.Dates(end)])

    % Display the legend and titles
    leg              = legend(Q,L) ;
    leg.Interpreter  = 'Latex' ;
    leg.FontSize     = FontSize ;
    title('Estimates of the COVID-19 reproduction number for optimal parameters $\lambda$','Interpreter','Latex')
    set(gca,'FontSize',FontSize,'ticklabelinterpreter','Latex')
    f5.Position      = [141 329 1033 314];
    if length(Q)     == 6
        leg.Position = [0.7103 0.2222 0.1948 0.6631] ;
    elseif length(Q) == 3
        leg.Position = [0.7103 0.5096 0.1948 0.3779] ;
    end
    
end