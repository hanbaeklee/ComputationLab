%% An RBC model with asset price, irreversibility, and endogenous labor supply with infinite Frisch elasticity
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
rbcassetirrendolaborinftyfrisch_ss;

% run the model with aggregate uncertainty
rbcassetirrendolaborinftyfrisch_bc;

% run the testers
rbcassetirrendolaborinftyfrisch_ee;
rbcassetirrendolaborinftyfrisch_monotonicity;

cd("../")
