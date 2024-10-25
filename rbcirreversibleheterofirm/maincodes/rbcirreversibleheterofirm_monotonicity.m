%% An RBC model with heterogeneous firms and irreversible investments
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
%load solutions
%=========================
ss = load('../solutions/rbcirreversibleheterofirm_ss.mat');
globalSol = load('../solutions/rbcirreversibleheterofirm_bc.mat');
load('../solutions/rbcirreversibleheterofirm_bc.mat');

%=========================
%backward solution
%=========================
iA = tsimpath;
vA  = vgridA(iA);
RHS = tjderiv;

%%

% the individual firm needs to be specified.
% the monotonicity result is robust over the choice of different individual
% firms.

iz = floor(pnumgridz/2);
ik = floor(pnumgridk/2);

tKsample = tK(burnin+1:pathlength-burnin);
RHSsample = RHS(:,:,burnin+1:pathlength-burnin);
tsimpathsample = tsimpath(burnin+1:pathlength-burnin);

for iAlocation = 1:pnumgridA
tempK = tKsample(tsimpathsample==iAlocation);
tempRHS = squeeze(RHSsample(ik,iz,tsimpathsample==iAlocation));
subplot(2,3,iAlocation);
scatter(tempK,tempRHS);
xlabel("K","FontSize",15);
ylabel("RHS of Euler","FontSize",15);
temptitle = append('A',num2str(iAlocation));
title(temptitle);
ax = gca;
ax.FontSize = 15; 
end
set(gcf, 'PaperPosition', [0 0 12 8]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [12 8]); %Set the paper to have width a and height b.
location = ['../figures/monotonicity.pdf'];
saveas(gcf, location);
