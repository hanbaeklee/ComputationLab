%% An RBC model with asset price, irreversibility, and endogenous labor supply
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
pFrisch         = 1.000;
pEeta           = 80.000;
pRiskAversion   = 2.000; 
pBbeta          = 0.960;

%=========================
% parameters - firms
%=========================
pAalpha         = 0.330;
pDdelta         = 0.100;
pPhi            = 0.975;

%%
%=========================
% steady-state equilibrium
%=========================
w           = (1-pAalpha)/(((1/pBbeta + pDdelta - 1))^(pAalpha/(1-pAalpha)));
K2L         = (w/(1-pAalpha))^(1/pAalpha);
L           = ((...
                w /( pEeta*( (((1-pAalpha)/w)^((1-pAalpha)/pAalpha) -pDdelta) *K2L)^pRiskAversion )...
               )^pFrisch)^(1/(1+pFrisch*pRiskAversion));

K           = K2L*L;
Y           = K^pAalpha*L^(1-pAalpha);
I           = pDdelta*K;
C           = Y - I;
J           = (C/(1-pBbeta));

%=========================  
% report
%=========================  
fprintf(' \n');

%=========================
% save
%=========================
dir = '../solutions/rbcassetirrendolabor_ss';
save(dir);
