% Load and display the synthetic piecewise linear reproduction coefficient
% used in Section IV of Pascal & Vaiter (2024) to generate synthetic
% nonstationary autoregressive Poisson time series.
%
% Synthetic ground truth reproduction coefficient designed to imitate the
% COVID-19 reproduction number evolution in France between October 8, 2021
% and August 3, 2022, corresponding to T = 300 days.
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

function [X, Y0] = synthetic_example(FontSize)

    % Input:   - FontSize: font size in the plots (optional, default FontSize = 22.5)
    %
    % Outputs: - X: piecewise linear reproduction coefficient used as ground truth in Section IV of Pascal & Vaiter (2024)
    %          - Y0: initial observation to generate the nonstationary autoregressive Poisson process

    if nargin < 1
        FontSize = 22.5 ;
    end
    
    %% LOAD SYNTHETIC REPRODUCTION COEFFICIENT

    load data/synthetic_example.mat
    

    %% DISPLAY THE REPRODUCTION COEFFICIENT WITH INDICATION OF Xt > 1

    % Color of the curve
    royal           = [0.2549019607843137, 0.4117647058823529, 0.8823529411764706] ;

    % Display
    f1 = figure(1); clf
    yline(1,'k-','linewidth',2)
    grid on
    hold on
    fill([1:111, fliplr(1:111)], [zeros(1,111), fliplr(2*ones(1,111))],royal,'FaceAlpha',0.15,'Edgecolor','none') ;
    fill([149:183, fliplr(149:183)], [zeros(1,35), fliplr(2*ones(1,35))],royal,'FaceAlpha',0.15,'Edgecolor','none') ;
    p1              = fill([236:283, fliplr(236:283)], [zeros(1,48), fliplr(2*ones(1,48))],royal,'FaceAlpha',0.15,'Edgecolor','none') ;
    p2              = plot(X,'-','linewidth',3,'color',royal) ;
    leg             = legend([p2,p1],'$\overline{\mathrm{X}}_t$','$\overline{\mathrm{X}}_t> 1$'); 
    leg.Interpreter = 'LaTex';
    leg.FontSize    = FontSize;
    leg.Position    = [0.1409 0.2144 0.1059 0.2847];
    xlabel('$t$','interpreter','latex')
    title('Piecewise linear reproduction coefficient with exponential growth periods','Interpreter','Latex')
    xlim([1,300])
    ylim([0,2])
    set(gca,'FontSize',FontSize,'ticklabelinterpreter','latex')
    f1.Position     = [141 329 1033 314];

end