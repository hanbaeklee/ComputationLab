%% An RBC model with asset price
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

%%
%=========================
% steady-state equilibrium
%=========================
K           = (pAalpha/(1/pBbeta + pDdelta - 1))^(1/(1-pAalpha));
Y           = K^pAalpha;
I           = pDdelta*K;
C           = Y - I;
J           = C/(1-pBbeta);

%=========================  
% report
%=========================  
fprintf(' \n');

%=========================
% save
%=========================
dir = '../solutions/rbcassetprice_ss';
save(dir);
