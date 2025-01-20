%% Basic NK model
%=========================    
% Unit period is a quarter
% Housekeeping
%=========================
clc;
clear variables;
close all; 
fnPath = '../functions';
addpath(fnPath);

%=========================
% households
%=========================
pFrisch         = 1.000;
pEta            = 2.050;
pRiskAversion   = 1.000; 
pBeta           = 0.990;

%=========================
% others
%=========================
pCalvo          = 0.660;
pElasticity     = 10.00;
pPsi            = (pElasticity - 1)*pCalvo/((1-pCalvo)*(1-pBeta*pCalvo));
pRhor           = 0.850; % 0.99;
pTaylorPi       = 1.500; % 1.50;
pTaylorY        = 0.125;
pPiTarget       = 0.02/4; % quarterly

%=========================
% auxiliary
%=========================
pNumGridA       = 1;
A               = 1;

%%
%=========================
% steady-state equilibrium
%=========================
error       = 10;
weightold   = 0.99;
tol_iter    = 1e-10;

b           = 0.00;
piAgg       = pPiTarget;
mc          = (1/pElasticity)*(pElasticity-1+pPsi*(1-pBeta)*(piAgg-pPiTarget)*(1+piAgg));
w           = mc*A;
Yf          = ((1/(pEta))*(pElasticity-1)/pElasticity)^(pFrisch/(1+pFrisch*pRiskAversion));

c           = (A*(w/pEta)^(pFrisch)*(1-(pPsi/2)*(piAgg-pPiTarget)^2))^(1/(1+pRiskAversion*pFrisch));
Y           = (1/(1-(pPsi/2)*(piAgg-pPiTarget)^2))*c;
n           = (w/(pEta*c^(pRiskAversion)))^(pFrisch);

%%
% additional allocations
% i           = (1+pPiTarget)*(1/pBeta)* ( (1+piAgg)/(1+pPiTarget) )^pTaylorPi-1;
i           = (1+pPiTarget)*(1/pBeta)* ( (1+piAgg)/(1+pPiTarget) )^pTaylorPi * (Y/Yf)^pTaylorY -1;

%=========================  
% report
%=========================  
fprintf(' \n');
fprintf('Consumption: %.10f \n', c);
fprintf('Marginal cost: %.10f \n', mc); 

%=========================
% save
%=========================
dir = '../solutions/nkrotem_ss';
save(dir);
