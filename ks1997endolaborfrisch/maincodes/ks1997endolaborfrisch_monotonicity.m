%% Krusell and Smith (1997) with endogenous labor supply
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
ss = load('../solutions/ks1997endolabor_ss.mat');
globalSol = load('../solutions/ks1997endolabor_bc.mat');
load('../solutions/ks1997endolabor_bc.mat');

%=========================
%backward solution
%=========================
iA = tsimpath;
vA  = vgridA(iA);

r = palpha.*vA.*(tK(1:end-1)./tsupplyL).^(palpha-1)-pdelta;
mr = zeros(size(mpolc));
for itrans = 1:pathlength
mr(:,:,itrans) = r(itrans);
end
RHS1 = (1+mr).*(1./mpolc);
RHS2 = (1./mpolc);


%%

% the individual household state needs to be specified.
% the monotonicity result is robust over the choice of different individual states.

iz = floor(pnumgridz/2);
ik = floor(pnumgridomega/2);

%RHS1
tKsample = tK(burnin+1:pathlength-burnin);
RHSsample = RHS1(:,:,burnin+1:pathlength-burnin);
tsimpathSample = tsimpath(burnin+1:pathlength-burnin);

figure;
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
set(gcf, 'PaperPosition', [0 0 9 4]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [9 4]); %Set the paper to have width a and height b.Grid off;
location = ['../figures/monotonicity1.pdf'];
saveas(gcf, location);


%RHS2
tKsample = tK(burnin+1:pathlength-burnin);
RHSsample = RHS2(:,:,burnin+1:pathlength-burnin);
tsimpathSample = tsimpath(burnin+1:pathlength-burnin);

figure;
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
set(gcf, 'PaperPosition', [0 0 9 4]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [9 4]); %Set the paper to have width a and height b.Grid off;
location = ['../figures/monotonicity2.pdf'];
saveas(gcf, location);