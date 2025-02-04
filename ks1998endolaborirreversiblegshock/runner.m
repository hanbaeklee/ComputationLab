%% Krusell and Smith (1998) with endogenous labor supply, investment irreversibility, and fiscal spending shock
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
ks1998endolaborirreversiblegshock_ss;

% run the model with aggregate uncertainty
ks1998endolaborirreversiblegshock_bc;

% run the testers
ks1998endolaborirreversiblegshock_ee;
ks1998endolaborirreversiblegshock_monotonicity;

cd("../")