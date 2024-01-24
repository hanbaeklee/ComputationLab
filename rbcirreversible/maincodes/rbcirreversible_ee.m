%% An RBC model with irreversible investment
% 2023.09.25
% Hanbaek Lee (hanbaeklee1@gmail.com)
% When you use the code, please cite the paper 
% "A Dynamically Consistent Global Nonlinear Solution 
% Method in the Sequence Space and Applications."
%=========================    
% this file is to compute the euler equation errors.
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
ss = load('../solutions/rbcirreversible_ss.mat');
globalSol = load('../solutions/rbcirreversible.mat');
load('../solutions/rbcirreversible.mat');

%%
%=========================
%backward solution
%=========================
iA = vSimPath;
iFuture = [(2:requiredTime)';requiredTime];
futureShock = vSimPath(iFuture);
vA  = vGridA(iA);
vKprime = [vK(2:end);vK(1)];

mTransException = zeros(size(vK));
for iTrans = 1:length(vK)
mTransException(iTrans,1) = mTransA(iA(iTrans),futureShock(iTrans));
end

tempV1 = 0;

for iAprime = 1:pNumGridA

    Aprime = vGridA(iAprime);

    candidate = vK(find(vSimPath==iAprime));
    candidateLocation = find(vSimPath==iAprime);
    candidate(candidateLocation>requiredTime-BURNIN) = [];
    candidate(candidateLocation<BURNIN) = [];
    candidateLocation(candidateLocation>requiredTime-BURNIN) = [];
    candidateLocation(candidateLocation<BURNIN) = [];
    [candidate,index] = sort(candidate);
    candidateLocation = candidateLocation(index);

    nLow = sum(repmat(candidate',length(vKprime),1)<vKprime,2);
    nLow(nLow<=1) = 1;
    nLow(nLow>=length(index)) = length(index)-1;
    nHigh = nLow+1;
    weightLow = (candidate(nHigh) - vKprime)./(candidate(nHigh)-candidate(nLow));
    weightLow(weightLow<0) = 0;
    weightLow(weightLow>1) = 1;
    
    rprime  = pAlpha.*Aprime.*vKprime.^(pAlpha-1) - pDelta;
    Lambdaprime = weightLow.*vLambda(candidateLocation(nLow )) + (1-weightLow).*vLambda(candidateLocation(nHigh ));

    tempV1 = tempV1 + (futureShock ~= iAprime).* pBeta.*...
                mTransA(iA,iAprime).*...
                (weightLow.*(1./vC(candidateLocation(nLow ))).^(pRiskAversion).*(1+rprime ) ...
           + (1-weightLow).*(1./vC(candidateLocation(nHigh))).^(pRiskAversion).*(1+rprime) ...
           - (1-pDelta)*Lambdaprime);

end

rfuture    = pAlpha*vGridA(futureShock).*vKprime.^(pAlpha-1) - pDelta;
tempV1     = tempV1 + pBeta*...
                mTransException.*((1./vC(iFuture)).^(pRiskAversion).*(1+rfuture) - (1-pDelta)*vLambda(iFuture)); 
tempV1     = tempV1 + vLambda;
tempC      = (1./tempV1).^(1/pRiskAversion);
EE         = (tempC - vC)./vC;

%%
ts1 = vK;%(vSimPath==3);
ts2 = EE;%(vSimPath==3);
% scatter(100*(log(ts1)-log(ss.K)),log10(abs(ts2)));
% xlabel("Aggregate capital K","FontSize",15);
figure;
scatter(1:length(ts2),log(abs(ts2)));
xlabel("Time (year)","FontSize",15);
ylabel("Absolute euler equation error (in log)","FontSize",15);
hold off;
location = ['../figures/ee_timeseries.pdf'];
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.Grid off;
saveas(gcf, location);

figure;
hist(log(abs(ts2)),100);
xlabel("Absolute euler equation error (in log)","FontSize",15);
ylabel("Distribution","FontSize",15);
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.Grid off;
location = ['../figures/ee_hist.pdf'];
saveas(gcf, location);

figure;
scatter(ts1,log(abs(ts2)));
xlabel("Aggregate capital stock","FontSize",15);
ylabel("Absolute euler equation error (in log)","FontSize",15);
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.Grid off;
location = ['../figures/ee.pdf'];
saveas(gcf, location);