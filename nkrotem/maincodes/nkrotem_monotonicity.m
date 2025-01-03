%% A canonical New Keynesian model (Rotemberg)
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
ss = load('../solutions/nkrotem_ss.mat');
globalSol = load('../solutions/nkrotem_bc.mat');
load('../solutions/nkrotem_bc.mat');

%=========================
%backward solution
%=========================
iA = vSimPath;
vA  = vGridA(iA);

RHS = (1./vc).^(pRiskAversion).*(1./(1+vPiAgg));
% for another inter-temporal optimality condition 
% RHS = (vc).^(-pRiskAversion) ...
%     .* (1+vPiAgg).*(vPiAgg-pPiTarget) ...
%     .* (vY);

%%

vSuffStatsample = vSuffStat(BURNIN+1:requiredTime-BURNIN);
RHSsample = RHS(BURNIN+1:requiredTime-BURNIN);
vsimpathSample = vSimPath(BURNIN+1:requiredTime-BURNIN);

for iAllocation = 1:5:pNumGridA
tempK = vSuffStatsample(vsimpathSample==iAllocation);
tempRHS = RHSsample(vsimpathSample==iAllocation);
subplot(3,5,(iAllocation+4)/5);
scatter(tempK,tempRHS);
xlabel("Suff stat","FontSize",15);
ylabel("Inter-temporal optimality","FontSize",15);
temptitle = append('A',num2str(iAllocation));
title(temptitle);
end
set(gcf, 'PaperPosition', [0 0 18 10]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [18 10]); %Set the paper to have width a and height b.Grid off;
location = ['../figures/monotonicity.pdf'];
saveas(gcf, location);
