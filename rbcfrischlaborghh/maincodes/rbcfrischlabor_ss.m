%% An RBC model with endogenous labor supply (Frisch elasticity-based)
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
pEeta           = 5.000;
pRiskAversion   = 1.000; 
pBbeta          = 0.990;

%=========================
% parameters - firms
%=========================
pAalpha         = 0.330;
pDdelta         = 0.025;

%=========================
% taxes
%=========================
% pTtauc          = 0.06;
% pTtauw          = 0.15;
% pTtaur          = 0.15;

%=========================
% aggregate shock
%=========================
% the steady state is with a trivial aggregate TFP
pNumGridA   = 1;
A           = 1;

%%
%=========================
% steady-state equilibrium
%=========================
r           = ((1/pBbeta) - 1);
K2L         = ((r+pDdelta)/(pAalpha*A))^(1/(pAalpha-1));
w           = K2L^pAalpha*A*(1-pAalpha);

L           = (w/pEeta)^pFrisch;
C           = (r*K2L + w)*L;
K           = K2L*L;
Y           = A*K^pAalpha*L^(1-pAalpha);


%=========================  
% report
%=========================  
fprintf(' \n');

%=========================
% save
%=========================
dir = '../solutions/rbcfrischlabor_ss';
save(dir);
