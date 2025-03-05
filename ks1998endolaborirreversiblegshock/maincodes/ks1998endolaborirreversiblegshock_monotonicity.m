%% Krusell and Smith (1998) with endogenous labor supply, investment irreversibility, and fiscal spending shock
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
ss = load('../solutions/ks1998endolaborirreversiblegshock_ss.mat');
globalSol = load('../solutions/ks1998endolaborirreversiblegshock_bc.mat');
load('../solutions/ks1998endolaborirreversiblegshock_bc.mat');

%=========================
%backward solution
%=========================
iA  = tsimpath;
vA  = vgridA(iA);
mr  = zeros(size(mpolc));
for itrans = 1:pathlength
  mr(:,:,itrans) = palpha*vgridA(tsimpath(itrans)).*(tK(itrans)./tsupplyL(itrans)).^(palpha-1)-pdelta;
end
RHS = ((1+mr)./mpolc - (1-pdelta)*mlambda);

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
subplot(3,7,iA);
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

%% 
% compute the spearman coefficients.
mtable = zeros(pnumgridz,pnumgrida,pnumgridA);
for iz = 1:pnumgridz
for ik = 1:pnumgrida

tKsample = tK(burnin+1:pathlength-burnin);
RHSsample = RHS(:,:,burnin+1:pathlength-burnin);
tsimpathSample = tsimpath(burnin+1:pathlength-burnin);

for iA = 1:pnumgridA
tempK = tKsample(tsimpathSample==iA);
tempRHS = squeeze(RHSsample(ik,iz,tsimpathSample==iA));

mtable(iz,ik,iA) = table2array(monotonicity( array2table([tempK,tempRHS]),"Method","rank"));

end

end
end

mean(mtable(:))
std(mtable(:))
min(mtable(:))