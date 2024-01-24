%% An RBC model with endogenous labor supply (Frisch elasticity-based)
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
ss = load('../solutions/rbcfrischlabor_ss.mat');
globalSol = load('../solutions/rbcfrischlabor_bc.mat');
load('../solutions/rbcfrischlabor_bc.mat');

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
    
    LLow    = ((1-pTtauw)*(1-pAalpha)*Aprime*vKprime.^pAalpha./(pEeta*(1+pTtauc)*vC(candidateLocation(nLow )).^pRiskAversion)).^(pFrisch/(1+pAalpha*pFrisch));
    LHigh   = ((1-pTtauw)*(1-pAalpha)*Aprime*vKprime.^pAalpha./(pEeta*(1+pTtauc)*vC(candidateLocation(nHigh )).^pRiskAversion)).^(pFrisch/(1+pAalpha*pFrisch));
    rLow    = pAalpha*Aprime*vKprime.^(pAalpha-1).*LLow.^(1-pAalpha) - pDdelta;
    rHigh   = pAalpha*Aprime*vKprime.^(pAalpha-1).*LHigh.^(1-pAalpha) - pDdelta;

    tempV1 = tempV1 + (futureShock ~= iAprime).* pBbeta .*...
                mTransA(iA,iAprime).*...
                (weightLow.*(1./vC(candidateLocation(nLow ))).^(pRiskAversion).*(1+(1-pTtaur).*rLow ) ...
           + (1-weightLow).*(1./vC(candidateLocation(nHigh))).^(pRiskAversion).*(1+(1-pTtaur).*rHigh) );

end
% for the realized future shock level on the simulated path
Lfuture    = ((1-pTtauw)*(1-pAalpha).*vGridA(futureShock).*vKprime.^pAalpha./(pEeta*(1+pTtauc)*vC(iFuture).^pRiskAversion)).^(pFrisch/(1+pAalpha*pFrisch));
rfuture    = pAalpha*vGridA(futureShock).*vKprime.^(pAalpha-1).*Lfuture.^(1-pAalpha) - pDdelta;
tempV1     = tempV1 + pBbeta*...
                mTransRealized.*(1./vC(iFuture)).^(pRiskAversion).*(1+(1-pTtaur).*rfuture); 

% update the allocations
tempC = ((1./tempV1).^(1/pRiskAversion));

EE         = (tempC - vC)./vC;

%%
ts1 = vK;%(vSimPath==3);
ts2 = EE;%(vSimPath==3);

figure;
scatter(1:length(ts2),log(abs(ts2)));
xlabel("Time (year)","FontSize",15);
ylabel("Euler equation error (in Log)","FontSize",15);
hold off;
location = ['../figures/ee_timeseries.pdf'];
saveas(gcf, location);

figure;
hist(log(abs(ts2)),100);
xlabel("Euler equation error (in Log)","FontSize",15);
ylabel("Distribution","FontSize",15);
location = ['../figures/ee_hist.pdf'];
saveas(gcf, location);

figure;
scatter(ts1,log(abs(ts2)));
xlabel("Aggregate capital stock","FontSize",15);
ylabel("Euler equation error (in Log)","FontSize",15);
location = ['../figures/ee.pdf'];
saveas(gcf, location);