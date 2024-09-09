clear all
close all
clc

addpath(genpath('/Users/bpascal/ownCloud/PhD_Barbara/Matlab/stein-piecewise-filtering'))
addpath(genpath(pwd))


%% GENERATE SYNTHETIC DATA ACCORDING TO SCALED POISSON (SP) MODEL

% Load and display the synthetic piecewise linear reproduction coefficient
[X, Y0]                 = synthetic_example() ;

% Model parameters
opts_SP.Y0              = Y0 ;  % initial observation
opts_SP.alpha           = 100; % Poisson scaling parameter (1 for standard Poisson)

% Sample a time series
[Y, Psi_Y, M]             = generate_synthetic_Poisson(X, opts_SP) ;
% - Y: synthetic observations
% - Psi_Y: memory functions evaluated in Y
% - M: parameters of the model
return
%% ESTIMATE THE REPRODUCTION COEFFICIENT THROUGH PIECEWISE LINEAR DENOISING 

% Parameter of the APURE unbiased risk estimate
alpha        = 100 ; % Poisson scaling parameter (default
N_MC         = 10 ; % number of Monte Carlo vectors of the robustified APURE estimates (default: 10)

% Setup for the grid search minimization of unbiased risk estimates
L            = 60 ; % number of explored values of the regularization parameter lambda
lambda_min   = 1e-2 ; % smallest lambda of the logarithmically spaced grid
lambda_max   = 1e4 ;% largest lambda of the logarithmically spaced grid

% Minimization of the APURE unbiased prediction risk estimate
opts_P.alpha      = alpha ;
opts_P.N_MC       = N_MC ;
opts_P.L          = L ;
opts_P.lambda_min = 1e-2 ;
opts_P.lambda_max = 1e4 ;
X_APURE_P         = APURE_Prediction(Y,Psi_Y,M,opts_P) ;

% Minimization of the APURE unbiased estimation risk estimate
opts_P.N_MC  = N_MC ;
opts_P.alpha = alpha ;
X_APURE_E    = APURE_Estimation(Y,Psi_Y,opts_E) ;