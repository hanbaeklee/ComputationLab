%% An RBC model with asset price and convex adjustment cost
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
pBbeta          = 0.960;

%=========================
% parameters - firms
%=========================
pAalpha         = 0.330;
pDdelta         = 0.100;
pMmu            = 0.800;

%%
%=========================
% steady-state equilibrium
%=========================
K           = (pAalpha/((1+pMmu*pDdelta)/pBbeta - (1-pDdelta) - pMmu/2 + (pMmu/2)*(1-pDdelta)^2))^(1/(1-pAalpha));
Y           = K^pAalpha;
I           = pDdelta*K;
C           = Y - I - (pMmu/2)*(pDdelta)^2*K;
J           = C/(1-pBbeta);
r           = pAalpha*K.^(pAalpha-1) + (1-pDdelta) + (pMmu/2) - (pMmu/2)*(1-pDdelta)^2;

%=========================  
% report
%=========================  
fprintf(' \n');

%=========================
% save
%=========================
dir = '../solutions/rbcassetcvxadjcost_ss';
save(dir);
