%% An RBC model with irreversible investment
% 2023.09.25
% Hanbaek Lee (hanbaeklee1@gmail.com)
% When you use the code, please cite the paper 
% "A Dynamically Consistent Global Nonlinear Solution 
% Method in the Sequence Space and Applications."
%=========================    
% this file is to compare the solutions across the different methods.
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
% load ss
%=========================
ss = load('../solutions/rbcirreversible_ss.mat');
globalSol = load('../solutions/rbcirreversible.mat');
occbinSol = load('../solutions/rbcirreversible_occbin.mat');
gdsgeSol  = load('../solutions/rbcirreversible_gdsge.mat');

globalSol_EE = load('../solutions/rbcirreversible.mat','EE');
occbinSol_EE = load('../solutions/rbcirreversible_occbin_EE.mat','EE');
occbinSollin_EE = load('../solutions/rbcirreversible_occbinlin_EE.mat','EE');
gdsgeSol_EE  = load('../solutions/rbcirreversible_gdsge_EE.mat','EE');

%%
%=========================
% variables
%=========================
%timing adjustment
requiredTime = globalSol.requiredTime;
occbinSol.vK = (1+occbinSol.k_p([1,1,1:requiredTime-2]))*ss.K;
occbinSol.vKlinear = (1+occbinSol.k_l([1,1,1:requiredTime-2]))*ss.K;
gdsgeSol.vK = gdsgeSol.SimuRslt.K([1,1:requiredTime-1])';

%%
%=========================
% plot
%=========================
figure;
% samplePeriod = 2000:3000;
samplePeriod = 4250:4300;
plot(samplePeriod,globalSol.vK(samplePeriod),"Color","blue","LineStyle","-","LineWidth",1.5); hold on;
plot(samplePeriod,occbinSol.vKlinear(samplePeriod),"Color","red","LineStyle","--","LineWidth",1.5); hold on;
plot(samplePeriod,occbinSol.vK(samplePeriod),"Color","black","LineStyle",":","LineWidth",1.5); hold on;
plot(samplePeriod,gdsgeSol.vK(samplePeriod),"Color","magenta","LineStyle","-.","LineWidth",1.5);
ylabel("Capital stock","FontSize",15);
xlabel("Time (year)","FontSize",15);
legend("rtm","linear","occbin","gdsge","fontsize",15,"location","northwest")
location = ['../figures/compsol.pdf'];
hold off;
box off
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.Grid off;
saveas(gcf, location);

%%
%=========================
% plot - EE
%=========================
figure;
% samplePeriod = 1000:4000;
samplePeriod = 4250:4300;
plot(samplePeriod,log(abs(globalSol_EE.EE(samplePeriod))),"Color","blue","LineStyle","-","LineWidth",1.5); hold on;
plot(samplePeriod,log(abs(occbinSollin_EE.EE(samplePeriod))),"Color","red","LineStyle","--","LineWidth",1.5); hold on;
plot(samplePeriod,log(abs(occbinSol_EE.EE(samplePeriod))),"Color","black","LineStyle",":","LineWidth",1.5); hold on;
plot(samplePeriod,log(abs(gdsgeSol_EE.EE(samplePeriod))),"Color","magenta","LineStyle","-.","LineWidth",1.5);
ylabel("Absolute euler equation error (in log)","FontSize",15);
xlabel("Time (year)","FontSize",15);
legend("rtm","linear","occbin","gdsge","fontsize",15,"location","southeast")
location = ['../figures/ee_tscomp.pdf'];
hold off;
box off;
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.Grid off;
saveas(gcf, location);
