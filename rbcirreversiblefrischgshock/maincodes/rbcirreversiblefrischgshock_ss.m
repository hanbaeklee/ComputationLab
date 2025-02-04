%% An RBC model with endogenous labor supply (Frisch elasticity-based), irreversible investment, and fiscal spending shock
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
pEeta           = 10.85;
pRiskAversion   = 1.000; 
pBbeta          = 0.960;
pPhi            = 0.975;%0.955;%0.500;%0.994;

%=========================
% parameters - firms
%=========================
pAalpha         = 0.330;
pDdelta         = 0.100;

%=========================
% taxes
%=========================
pOmega          = 0.200;

%=========================
% aggregate shock
%=========================
% the steady state is with a trivial aggregate TFP
pNumGridA       = 1;
A               = 1;

%%
%=========================
% steady-state equilibrium
%=========================
error       = 0.2;
weightOld   = 0.8;

C           = 0.7;
while error>1e-10
r           = (1/pBbeta) - 1;
K2L         = ((r+pDdelta)/(pAalpha*A))^(1/(pAalpha-1));
w           = K2L^pAalpha*A*(1-pAalpha);
L           = (w/(pEeta*C^pRiskAversion))^pFrisch;
K           = K2L*L;
I           = K*pDdelta;
Y           = A*K^pAalpha*L^(1-pAalpha);
T           = Y*pOmega;
G           = T;
Lambda      = 0;
Cnew        = Y - I - T; 
error       = abs(Cnew - C);
C           = weightOld*C + (1-weightOld)*Cnew;
end

%=========================  
% report
%=========================  
fprintf(' \n');

%=========================
% save
%=========================
dir = '../solutions/rbcirreversiblefrischgshock_ss';
save(dir);
