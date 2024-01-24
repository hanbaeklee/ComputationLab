%% An RBC model with heterogeneous firms and irreversible investments
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
rbcirreversibleheterofirm_ss;

% run the model with aggregate uncertainty
rbcirreversibleheterofirm_bc;

% run the testers
rbcirreversibleheterofirm_monotonicity;

% comparison with the representative-agent model
rbcirreversiblecomp_heterorep

cd("../")
