%% Krusell and Smith (1998) with endogenous labor supply and convex adjustment cost
% 2023.09.25
% Hanbaek Lee (hanbaeklee1@gmail.com)
% When you use the code, please cite the paper 
% "A Dynamically Consistent Global Nonlinear Solution 
% Method in the Sequence Space and Applications."
%=========================    
% this file is to check the monotonicity of value/policy function's
% along with the endogenous state in equilibrium.
%=========================    
%=========================    
% housekeeping
%=========================
% clc;
% clear variables;
% close all; 
% fnPath = './functions';
% addpath(fnPath);

%=========================
%load ss
%=========================
ss = load('../solutions/ks1998endolaborfrischadjcost_ss.mat');
globalSol = load('../solutions/ks1998endolaborfrischadjcost_bc.mat');
load('../solutions/ks1998endolaborfrischadjcost_bc.mat');

%=========================
%backward solution
%=========================
iA = tsimpath;
vA  = vgridA(iA);
RHS = mpolc;

%%

% the individual household state needs to be specified.
% the monotonicity result is robust over the choice of different individual states.

iz = floor(pnumgridz/2);
ik = floor(pnumgrida/2);

tKsample = tK(burnin+1:pathlength-burnin);
RHSsample = RHS(:,:,burnin+1:pathlength-burnin);
tsimpathSample = tsimpath(burnin+1:pathlength-burnin);

for iA = 1:pnumgridA
tempK = tKsample(tsimpathSample==iA);
tempRHS = squeeze(RHSsample(ik,iz,tsimpathSample==iA));
subplot(1,2,iA);
scatter(tempK,tempRHS);
xlabel("K","FontSize",15);
ylabel("RHS of Euler","FontSize",15);
temptitle = append('A',num2str(iA));
title(temptitle);
end
set(gcf, 'PaperPosition', [0 0 18 10]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [18 10]); %Set the paper to have width a and height b.Grid off;
location = ['../figures/monotonicity.pdf'];
saveas(gcf, location);
