%% An RBC model with endogenous labor supply (Frisch elasticity-based)
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
rbcfrischlabor_ss;

% run the model with aggregate uncertainty
rbcfrischlabor_bc;

% run the testers
rbcfrischlabor_ee;
rbcfrischlabor_monotonicity;

cd("../")

