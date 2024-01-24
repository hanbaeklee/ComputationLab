%% A canonical RBC model
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
pTtheta         = 1.7517;
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
pTtauc          = 0.00;%0.06;
pTtauw          = 0.00;%0.15;
pTtaur          = 0.00;%0.15;

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
r           = ((1/pBbeta) - 1)*(1/(1-pTtaur));
K2L         = ((r+pDdelta)/(pAalpha*A))^(1/(pAalpha-1));
w           = K2L^pAalpha*A*(1-pAalpha);

L           = ((1-pTtauw)*w)/(pTtheta*(1+pTtauc)*(r*K2L+w)+(1-pTtauw)*w);
C           = (r*K2L + w)*L;
K           = K2L*L;
Y           = A*K^pAalpha*L^(1-pAalpha);
T           = pTtaur*r*K + pTtauw*w*L + pTtauc*C;

%=========================  
% report
%=========================  
fprintf(' \n');

%=========================
% save
%=========================
dir = '../solutions/rbc_ss';
save(dir);
