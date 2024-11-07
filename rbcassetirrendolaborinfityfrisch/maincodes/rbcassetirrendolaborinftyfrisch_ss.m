%% An RBC model with asset price, irreversibility, and endogenous labor supply with infinite Frisch elasticity
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
pbeta           = 0.977;
peta            = 2.400;

%=========================
% parameters - firms
%=========================
palpha          = 0.256;
pgamma          = 0.640;
pdelta          = 0.069;
pphi            = 0.975;

%=========================
% numerical parameters
%=========================
tol_ge          = 1e-10;
weightold       = 0.800;

%=========================
% initial guess
%=========================
eq.p               = peta;
eq.k               = 0.3;

%%
%=========================
% steady-state equilibrium
%=========================
error2 = 10;
while error2>tol_ge

eq.w            = peta/eq.p;
eq.l            = (pgamma*eq.k^palpha/eq.w)^(1/(1-pgamma));
eq.i            = pdelta*eq.k;
eq.y            = eq.k^palpha*eq.l^pgamma;
eq.c            = eq.k^palpha*eq.l^pgamma - eq.i;
eq.j            = (eq.y - eq.w*eq.l - eq.i)/(1-pbeta);
eq.r            = (1-pgamma)*(pgamma/eq.w).^(pgamma/(1-pgamma))...
                .* (palpha/(1-pgamma)).*eq.k.^(palpha/(1-pgamma)-1) ...
                + (1-pdelta);

pnew            = 1/eq.c;
knew            = ( (1/pbeta - (1-pdelta))...
                  /((1-pgamma)*(pgamma/eq.w).^(pgamma/(1-pgamma))...
                    .* (palpha/(1-pgamma))) )^(1/(palpha/(1-pgamma)-1)) ;
                 
error2          = max(abs([eq.p-pnew,eq.k-knew]));
eq.p            = weightold * eq.p + (1-weightold) * pnew;
eq.k            = weightold * eq.k + (1-weightold) * knew;


end

%=========================  
% report
%=========================  
fprintf(' \n');

%=========================
% save
%=========================
dir = '../solutions/rbcassetirrendolaborinftyfrisch_ss';
save(dir);
