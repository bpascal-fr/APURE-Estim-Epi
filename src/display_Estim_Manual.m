% Display the estimates of the reproduction coefficient either 
% corresponding to the Maximum Likelihood Estimate under the autoregressive
% Poisson model and to the Penalized Kullback-Leibler variational estimate
% favoring piecewise linearity.
%
% B. Pascal
% April, 2025

function display_Estim_Manual(X, FontSize)

    % Inputs:  - X: estimated piecewise linear reproduction coefficients
    %            - GT: ground truh (if available)
    %            - ML: estimated maximum likelihood reproduction coefficient
    %            - PKL: estimated penalized Kullback-Leibler reproduction coefficient
    %            - Dates: abstract dates in datetime format or time indices for display (optional, by default 1 to T)
    %          - FontSize: font size used in the plots (optional, by default 22.5)
    %

    if nargin < 2
        FontSize = 22.5 ;
    end

    if ~isfield(X,'Dates');    X.Dates = 1:length(X.ML);    end

    % Customized colors
    royal  = [0.2549019607843137, 0.4117647058823529, 0.8823529411764706] ;  
    cloud  = [0.65,0.65,0.65] ;
    brown  = [137,81,41]/255 ;
    
    % Display the estimates and ground truth if available
    f3  = figure(3); clf

    % Reference vector full of ones
    plot(X.Dates,ones(size(X.ML)),'-k')
    hold on
    grid on
    Q    = [] ;
    iEst = 1 ; 
    M    = 1 ;

    % Display ground truth if provided
    if isfield(X,'GT')
        q            = plot(X.Dates,X.GT,'linewidth',3,'color',royal) ;
        Q            = [Q, q] ;
        L{iEst}      = '$\overline{\mathrm{X}}_t$' ;
        M            = max(M, max(X.GT)) ;
        iEst         = iEst + 1 ;
    end

    % Display Maximum Likelihood Estimate
    if isfield(X,'ML')
        q       = plot(X.Dates,X.ML,'linewidth',1,'color',cloud) ;
        Q       = [Q, q] ;
        L{iEst} = '$\mathrm{X}_t^{\mathrm{ML}}$' ;
        M       = max(M, max(X.ML)) ;
        iEst    = iEst + 1 ;
    end

    % Display the estimate obtained minimizing the true risk if available
    if isfield(X,'PKL')
        q            = plot(X.Dates,X.PKL,'linewidth',3,'color',brown) ;
        Q            = [Q, q] ;
        L{iEst}      = '$\widehat{\mathrm{X}}_t(\mathrm{Y};\lambda)$' ;
        M            = max(M, max(X.PKL)) ;
        iEst         = iEst + 1 ;
    end


    % Limits of the y-axis
    ylim([0, 1.1 * M])

    % Display the legend and titles
    l3               = legend(Q,L) ;
    l3.Interpreter   = 'Latex' ;
    if length(Q) == 3
        l3.Position  = [0.7868 0.5986 0.1183 0.2868] ;
    elseif length(Q) == 2
        l3.Position  =  [0.8245 0.6959 0.0806 0.1895] ;
    end
    l3.FontSize      = FontSize ;
    l3.Color         = 'None' ;
    title('Estimates of the reproduction coefficients with no or manual tuning','Interpreter','Latex')
    set(gca,'FontSize',FontSize,'ticklabelinterpreter','Latex','Color','none')
    f3.Position      = [141 329 1033 314];

end