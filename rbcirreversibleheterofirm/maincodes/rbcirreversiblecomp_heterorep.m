%% An RBC model with heterogeneous firms and irreversible investments
% 2023.09.25
% Hanbaek Lee (hanbaeklee1@gmail.com)
% When you use the code, please cite the paper 
% "A Dynamically Consistent Global Nonlinear Solution 
% Method in the Sequence Space and Applications."
%=========================    
% this file is to compare the heterogeneous-agent model with the
% representative-agent model
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
ss      = load('../solutions/rbcirreversibleheterofirm_ss.mat');
hetero  = load('../solutions/rbcirreversibleheterofirm_bc.mat');
rep     = load('../../rbcassetirrendolaborinfityfrisch/solutions/rbcassetirrendolaborinftyfrisch_bc.mat');
heterofrless  = load('../solutions/rbcirreversibleheterofirm_bc_frictionless.mat');
repfrless     = load('../../rbcassetirrendolaborinfityfrisch/solutions/rbcassetirrendolaborinftyfrisch_bc_frictionless.mat');

sampleperiod = 600:1000;
burnin = hetero.burnin;

%%
%=========================
%compare capital dynamics
%=========================
% define the series
ts1     = hetero.tK(sampleperiod); 
ts2     = rep.tk(sampleperiod); 

% normalization
ts1     = log(ts1) - log(mean(ts1));
ts2     = log(ts2) - log(mean(ts2));

yRange = [min([ts1;ts2])*1.1,max([ts1;ts2])*1.1];
xRange = [min(sampleperiod-hetero.burnin),max(sampleperiod-hetero.burnin)];

figure;
plot(sampleperiod-hetero.burnin,ts1,'linewidth',3);hold on;
plot(sampleperiod-hetero.burnin,ts2,'-.','linewidth',2,'Color',"red");
ylim(yRange);
xlim(xRange)
xlabel('Time (year)',"FontSize",15);
legend("Hetero","Rep.","location","best","FontSize",15);
hold off;
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.Grid off;

% save
location = ['../figures/comp_heterorep_capital','.pdf'];
saveas(gcf, location);


%%
%=========================
%compare capital dynamics - frictionless
%=========================
% define the series
ts1     = heterofrless.tK(sampleperiod); 
ts2     = repfrless.tk(sampleperiod); 
ts3     = heterofrless.recovered(sampleperiod-burnin)';
% normalization
ts1     = log(ts1) - log(mean(ts1));
ts2     = log(ts2) - log(mean(ts2));
ts3     = log(ts3) - log(mean(ts3));

yRange = [min([ts1;ts2;ts3])*1.1,max([ts1;ts2;ts3])*1.1];
xRange = [min(sampleperiod-hetero.burnin),max(sampleperiod-hetero.burnin)];

figure;
plot(sampleperiod-hetero.burnin,ts1,'linewidth',3);hold on;
plot(sampleperiod-hetero.burnin,ts2,'-.','linewidth',2,'Color',"red");
plot(sampleperiod-hetero.burnin,ts3,'--','linewidth',2,'Color',"black");
ylim(yRange);
xlim(xRange)
xlabel('Time (year)',"FontSize",15);
legend("Hetero","Rep.","Linear LoM","location","best","FontSize",15);
hold off;
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.Grid off;

% save
location = ['../figures/comp_heterorep_capital_frictionless','.pdf'];
saveas(gcf, location);

%%
%=========================
%compare consumption dynamics
%=========================
% define the series
ts1     = hetero.eq.c(sampleperiod); 
ts2     = rep.tc(sampleperiod); 

% normalization
ts1     = log(ts1) - log(mean(ts1));
ts2     = log(ts2) - log(mean(ts2));

yRange = [min([ts1;ts2])*1.1,max([ts1;ts2])*1.1];
xRange = [min(sampleperiod-hetero.burnin),max(sampleperiod-hetero.burnin)];

figure;
plot(sampleperiod-hetero.burnin,ts1,'linewidth',3);hold on;
plot(sampleperiod-hetero.burnin,ts2,'-.','linewidth',2,'Color',"red");
ylim(yRange);
xlim(xRange)
xlabel('Time (year)',"FontSize",15);
legend("Hetero","Rep.","location","best","FontSize",15);
hold off;
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.Grid off;

% save
location = ['../figures/comp_heterorep_consumption','.pdf'];
saveas(gcf, location);

%%
%=========================
%compare consumption dynamics - frictionless
%=========================
% define the series
ts1     = heterofrless.eq.c(sampleperiod); 
ts2     = repfrless.tc(sampleperiod); 

% normalization
ts1     = log(ts1) - log(mean(ts1));
ts2     = log(ts2) - log(mean(ts2));

yRange = [min([ts1;ts2])*1.1,max([ts1;ts2])*1.1];
xRange = [min(sampleperiod-hetero.burnin),max(sampleperiod-hetero.burnin)];

figure;
plot(sampleperiod-hetero.burnin,ts1,'linewidth',3);hold on;
plot(sampleperiod-hetero.burnin,ts2,'-.','linewidth',2,'Color',"red");
ylim(yRange);
xlim(xRange)
xlabel('Time (year)',"FontSize",15);
legend("Hetero","Rep.","location","best","FontSize",15);
hold off;
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.Grid off;

% save
location = ['../figures/comp_heterorep_consumption_frictionless','.pdf'];
saveas(gcf, location);

%%
%=========================
%compare output dynamics
%=========================
% define the series
ts1     = hetero.eq.y(sampleperiod)'; 
ts2     = rep.ty(sampleperiod); 

% normalization
ts1     = log(ts1) - log(mean(ts1));
ts2     = log(ts2) - log(mean(ts2));

yRange = [min([ts1;ts2])*1.1,max([ts1;ts2])*1.1];
xRange = [min(sampleperiod-hetero.burnin),max(sampleperiod-hetero.burnin)];

figure;
plot(sampleperiod-hetero.burnin,ts1,'linewidth',3);hold on;
plot(sampleperiod-hetero.burnin,ts2,'-.','linewidth',2,'Color',"red");
ylim(yRange);
xlim(xRange)
xlabel('Time (year)',"FontSize",15);
legend("Hetero","Rep.","location","best","FontSize",15);
hold off;
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.Grid off;

% save
location = ['../figures/comp_heterorep_output','.pdf'];
saveas(gcf, location);


%%
%=========================
%compare output dynamics - frictionless
%=========================
% define the series
ts1     = heterofrless.eq.y(sampleperiod)'; 
ts2     = repfrless.ty(sampleperiod); 

% normalization
ts1     = log(ts1) - log(mean(ts1));
ts2     = log(ts2) - log(mean(ts2));

yRange = [min([ts1;ts2])*1.1,max([ts1;ts2])*1.1];
xRange = [min(sampleperiod-hetero.burnin),max(sampleperiod-hetero.burnin)];

figure;
plot(sampleperiod-hetero.burnin,ts1,'linewidth',3);hold on;
plot(sampleperiod-hetero.burnin,ts2,'-.','linewidth',2,'Color',"red");
ylim(yRange);
xlim(xRange)
xlabel('Time (year)',"FontSize",15);
legend("Hetero","Rep.","location","best","FontSize",15);
hold off;
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.Grid off;

% save
location = ['../figures/comp_heterorep_output_frictionless','.pdf'];
saveas(gcf, location);

%%
%=========================
%compare investment dynamics
%=========================
% define the series
ts1     = hetero.eq.i(sampleperiod)'; 
ts2     = rep.ti(sampleperiod); 

% normalization
ts1     = log(ts1) - log(mean(ts1));
ts2     = log(ts2) - log(mean(ts2));

yRange = [min([ts1;ts2])*1.1,max([ts1;ts2])*1.1];
xRange = [min(sampleperiod-hetero.burnin),max(sampleperiod-hetero.burnin)];

figure;
plot(sampleperiod-hetero.burnin,ts1,'linewidth',3);hold on;
plot(sampleperiod-hetero.burnin,ts2,'-.','linewidth',2,'Color',"red");
ylim(yRange);
xlim(xRange)
xlabel('Time (year)',"FontSize",15);
legend("Hetero","Rep.","location","best","FontSize",15);
hold off;
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.Grid off;

% save
location = ['../figures/comp_heterorep_investment','.pdf'];
saveas(gcf, location);


%%
%=========================
%compare investment dynamics - frictionless
%=========================
% define the series
ts1     = heterofrless.eq.i(sampleperiod)'; 
ts2     = repfrless.ti(sampleperiod); 

% normalization
ts1     = log(ts1) - log(mean(ts1));
ts2     = log(ts2) - log(mean(ts2));

yRange = [min([ts1;ts2])*1.1,max([ts1;ts2])*1.1];
xRange = [min(sampleperiod-hetero.burnin),max(sampleperiod-hetero.burnin)];

figure;
plot(sampleperiod-hetero.burnin,ts1,'linewidth',3);hold on;
plot(sampleperiod-hetero.burnin,ts2,'-.','linewidth',2,'Color',"red");
ylim(yRange);
xlim(xRange)
xlabel('Time (year)',"FontSize",15);
legend("Hetero","Rep.","location","best","FontSize",15);
hold off;
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.Grid off;

% save
location = ['../figures/comp_heterorep_investment_frictionless','.pdf'];
saveas(gcf, location);


%%
%=========================
%compare business cycle statistics
%=========================
fprintf('Comparison for the log series');
fprintf(' \n');
disp(['volatility - y: ',num2str(round(std(log(hetero.eq.y)),3)),' vs ',num2str(round(std(log(rep.ty)),3))]);
fprintf(' \n');
disp(['volatility - c: ',num2str(round(std(log(hetero.eq.c)),3)),' vs ',num2str(round(std(log(rep.tc)),3))]);
fprintf(' \n');
disp(['volatility - i: ',num2str(round(std(log(hetero.eq.i)),3)),' vs ',num2str(round(std(log(rep.ti)),3))]);
fprintf(' \n');
fprintf(' \n');
disp(['skewness - y: ',num2str(round(skewness(log(hetero.eq.y)),3)),' vs ',num2str(round(skewness(log(rep.ty)),3))]);
fprintf(' \n');
disp(['skewness - c: ',num2str(round(skewness(log(hetero.eq.c)),3)),' vs ',num2str(round(skewness(log(rep.tc)),3))]);
fprintf(' \n');
disp(['skewness - i: ',num2str(round(skewness(log(hetero.eq.i)),3)),' vs ',num2str(round(skewness(log(rep.ti)),3))]);


