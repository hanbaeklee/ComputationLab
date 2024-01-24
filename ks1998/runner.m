%% Solving the model in Krusell and Smith (1998)
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
ks1998_ss;

% run the model with aggregate uncertainty
ks1998_bc;

% run the testers
ks1998_ee;
ks1998_monotonicity;

cd("../")