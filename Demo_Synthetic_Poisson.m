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
[Y, Psi, M]             = generate_synthetic_Poisson(X, opts_SP) ;
% - Y: synthetic observations
% - Psi: memory functions evaluated in Y
% - M: parameters of the model
