clear all
close all
clc

addpath(genpath('/Users/bpascal/ownCloud/PhD_Barbara/Matlab/stein-piecewise-filtering'))
addpath(genpath(pwd))

%% LOAD AND DISPLAY NEW INFECTION COUNTS

% List all countries referenced in JHU repository
list             = 1;
AllCountries(list)
% list = 0: do not display countries
% list = 1: display all countries names

% Countries to monitor
User_Country                = "France" ; 

% Time period
opts_load.LastDay           = '2023-02-05'; % Last day of the time period in format 'YYYY-MM-DD'. If not provided set to March 9, 2023.
opts_load.W                 = 70;           % Length of the time period in weeks. If -1 or not provided entire available time period.

% Load data
[Z, Phi_Z, M]               = load_JHU_Weekly(User_Country,opts_load);

%% MAXIMUM LIKELIHOOD ESTIMATE

X_MLE             = X_MaxLikelihood(Z,Phi_Z) ; 

%% ESTIMATE THE REPRODUCTION COEFFICIENT THROUGH PIECEWISE LINEAR DENOISING 

% Parameter of the APURE unbiased risk estimate
N                 = 1 ;   % number of Monte Carlo vectors of the robustified APURE estimates (default: 10)

% Setup for the grid search minimization of unbiased risk estimates
L                 = 2 ;   % number of explored values of the regularization parameter lambda
lambda_min        = 1e-2 ; % smallest lambda of the logarithmically spaced grid
lambda_max        = 1e4 ;  % largest lambda of the logarithmically spaced grid

% Minimization of the APURE unbiased prediction risk estimate
opts.N          = N ;
opts.L          = L ;
opts.lambda_min = 1e-2 ;
opts.lambda_max = 1e4 ;

% Minimization of the APURE unbiased prediction risk estimate
[R_P, lambda_P, oracle_P] = APURE_Prediction(Z,Phi_Z,M,opts) ;

% Minimization of the APURE unbiased estimation risk estimate
[R_E, lambda_E, oracle_E] = APURE_Estimation(Z,Phi_Z,M,opts) ;

%% DISPLAY ORACLES AND ESTIMATED REPRODUCTION COEFFICIENTS

% Minimizing oracles of the prediction quadratic risk
display_Covid_Prediction(R_P, lambda_P, oracle_P)

% Minimizing oracles of the estimation quadratic risk
display_Covid_Estimation(R_E, lambda_E, oracle_E)

% Compare the minimization of the prediction and estimation risks
compare_Covid(R_E, R_P)

