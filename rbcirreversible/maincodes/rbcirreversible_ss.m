%% An RBC model with irreversible investment
% 2023.09.25
% Hanbaek Lee (hanbaeklee1@gmail.com)
% When you use the code, please cite the paper 
% "A Dynamically Consistent Global Nonlinear Solution 
% Method in the Sequence Space and Applications."
%=========================    
% this file is to compute the steady-state equilibrium.
%=========================    
%=========================    
% housekeeping
%=========================
clc;
clear variables;
close all; 
fnPath = '../functions';
addpath(fnPath);

%=========================
% parameters - households
%=========================
pRiskAversion   = 2.000; 
pBeta          = 0.960;

%=========================
% parameters - firms
%=========================
pAlpha         = 0.330;
pDelta         = 0.100;
pPhi            = 0.975;

%%
%=========================
% steady-state equilibrium
%=========================
error       = 10;
K           = (pAlpha/(1/pBeta + pDelta - 1))^(1/(1-pAlpha));
Y           = K^pAlpha;
I           = pDelta*K;
C           = Y - I;

%=========================  
% report
%=========================  
fprintf(' \n');

%=========================
% save
%=========================
dir = '../solutions/rbcirreversible_ss';
save(dir);
