%% An RBC model with asset price and convex adjustment cost under TFP and adjustment friction fluctuations
% 2023.09.25
% Hanbaek Lee (hanbaeklee1@gmail.com)
% When you use the code, please cite the paper 
% "A Dynamically Consistent Global Nonlinear Solution 
% Method in the Sequence Space and Applications."
%=========================    
% this file is to run the whole codes.
%=========================    
%=========================
% housekeeping
%=========================
clc;
clear variables;
close all;
fnpath = './functions';
addpath(fnpath);

%%
%=========================
% run the codes
%=========================
cd("./maincodes/")

% obtain the steady state
rbcassetcvx2shocks_ss;

% run the model with aggregate uncertainty
rbcassetcvx2shocks_bc;

% run the testers
rbcassetcvx2shocks_ee;
rbcassetcvx2shocks_monotonicity;

cd("../")