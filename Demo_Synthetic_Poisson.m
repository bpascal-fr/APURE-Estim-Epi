clear all
close all
clc

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

X_ML              = X_MaxLikelihood(Y, Psi_Y) ; 


%% PENALIZED KULLBACK-LEIBLER ESTIMATE WITH MANUAL CHOICE OF THE REGULARIZATION PARAMETER

% Manual choice of the regularization parameter
lambda            = 50 ;

% Minimization of the penalized Kullback-Leibler functional
X_PKL             = X_Penalized(Y, Psi_Y, lambda) ;


%% DISPLAY THE MAXIMUM LIKELIHOOD AND VARIATIONAL ESTIMATE WITH MANUALLY TUNED REGULARIZATION PARAMETER

% Store the quantities to be displayed
X_Manual.GT       = X ;
X_Manual.ML       = X_ML ;
X_Manual.PKL      = X_PKL ;

% Plot the ground truth, if available, and manual estimates
display_Estim_Manual(X_Manual)

%% ESTIMATE THE REPRODUCTION COEFFICIENT THROUGH PIECEWISE LINEAR DENOISING 

% Minimization of the APURE unbiased prediction risk estimate
[X_P, lambda_P, oracle_P] = APURE_Prediction(Y,Psi_Y,M) ;

% Minimization of the APURE unbiased estimation risk estimate
[X_E, lambda_E, oracle_E] = APURE_Estimation(Y,Psi_Y,M) ;


%% DISPLAY ORACLES AND ESTIMATED REPRODUCTION COEFFICIENTS

% Minimizing oracles of the prediction quadratic risk
display_APURE_Prediction(X_P, lambda_P, oracle_P)

% Minimizing oracles of the estimation quadratic risk
display_APURE_Estimation(X_E, lambda_E, oracle_E)

% Compare the minimization of the prediction and estimation risks
compare_APURE(X_E, X_P)

