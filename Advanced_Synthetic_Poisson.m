clear all
close all
clc

addpath(genpath('/Users/bpascal/ownCloud/PhD_Barbara/Matlab/stein-piecewise-filtering'))
addpath(genpath(pwd))


%% GENERATE SYNTHETIC DATA ACCORDING TO SCALED POISSON MODEL

% Load and display the synthetic piecewise linear reproduction coefficient
[X, Y0]                 = synthetic_example() ;

% Model parameters
sett.Y0                 = Y0 ;  % initial observation
sett.alpha              = 1000; % Poisson scaling parameter (1 for standard Poisson)

% Sample a time series
[Y, Psi_Y, M]           = generate_synthetic_Poisson(X, sett) ;
% - Y: synthetic observations
% - Psi_Y: memory functions evaluated in Y
% - M: parameters of the model

%% MAXIMUM LIKELIHOOD ESTIMATE

X_MLE             = X_MaxLikelihood(Y, Psi_Y) ; 

%% ESTIMATE THE REPRODUCTION COEFFICIENT THROUGH PIECEWISE LINEAR DENOISING 

% Parameter of the APURE unbiased risk estimate
alpha             = 1000 ; % Poisson scaling parameter (default: 1 for standard Poisson)
N                 = 5  ;   % number of Monte Carlo vectors of the robustified APURE estimates (default: 10)

% Setup for the grid search minimization of unbiased risk estimates
L                 = 30 ;   % number of explored values of the regularization parameter lambda (default: 60)
lambda_min        = 1e-2 ; % smallest lambda of the logarithmically spaced grid
lambda_max        = 1e4 ;  % largest lambda of the logarithmically spaced grid

% Minimization of the APURE unbiased prediction risk estimate
opts_P.alpha      = alpha ;
opts_P.N          = N ;
opts_P.L          = L ;
opts_P.lambda_min = 1e-2 ;
opts_P.lambda_max = 1e4 ;
[X_P, lambda_P, oracle_P] = APURE_Prediction(Y,Psi_Y,M,opts_P) ;
% - X_P.GT: ground truth (if provided)
% - X_P.MLE: estimated maximum likelihood reproduction coefficient
% - X_P.RISK: estimated reproduction coefficient with regularization parameter minimizing true prediction risk (if ground truth available)
% - X_P.APURE: estimated reproduction coefficient with regularization parameter minimizing APURE prediction risk estimate
% - X_P.Dates: abstract dates in datetime format for display

% Minimization of the APURE unbiased estimation risk estimate
opts_E.alpha      = alpha ;
opts_E.N          = N ;
opts_E.L          = L ;
opts_E.lambda_min = 1e-2 ;
opts_E.lambda_max = 1e4 ;
[X_E, lambda_E, oracle_E] = APURE_Estimation(Y,Psi_Y,M,opts_E) ;
% - X_E.GT: ground truth (if provided)
% - X_E.MLE: estimated maximum likelihood reproduction coefficient
% - X_E.RISK: estimated reproduction coefficient with regularization parameter minimizing true estimation risk (if ground truth available)
% - X_E.APURE: estimated reproduction coefficient with regularization parameter minimizing APURE estimation risk estimate
% - X_E.Dates: abstract dates in datetime format for display

%% DISPLAY ORACLES AND ESTIMATED REPRODUCTION COEFFICIENTS

% Minimizing oracles of the prediction quadratic risk
display_APURE_Prediction(X_P, lambda_P, oracle_P)

% Minimizing oracles of the estimation quadratic risk
display_APURE_Estimation(X_E, lambda_E, oracle_E)

% Compare the minimization of the prediction and estimation risks
compare_APURE(X_E, X_P)

