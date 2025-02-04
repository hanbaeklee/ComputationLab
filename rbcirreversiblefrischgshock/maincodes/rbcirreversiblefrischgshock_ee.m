%% An RBC model with endogenous labor supply (Frisch elasticity-based), irreversible investment, and fiscal spending shock
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
ss = load('../solutions/rbcirreversiblefrischgshock_ss.mat');
globalSol = load('../solutions/rbcirreversiblefrischgshock_bc.mat');
load('../solutions/rbcirreversiblefrischgshock_bc.mat');

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

% declare an empty object tempV1 that will carry the cumulatively summed expected
% values.
tempV1 = 0;
for iAprime = 1:pNumGridA

    Aprime = vGridA(iAprime);

    % find a period where the future shock realization is the same as
    % iAprime and the capital stock is closest to vKprime from below and above.
    candidate = vK(find(vSimPath==iAprime)); % iso-shock periods
    candidateLocation = find(vSimPath==iAprime); % iso-shock period locations
    candidate(candidateLocation>requiredTime-BURNIN) = []; % last burnin periods cannot be a candidate
    candidate(candidateLocation<BURNIN) = [];  % initial burnin periods cannot be a candidate
    candidateLocation(candidateLocation>requiredTime-BURNIN) = []; % last burnin periods cannot be a candidate
    candidateLocation(candidateLocation<BURNIN) = [];  % initial burnin periods cannot be a candidate
    [candidate,index] = sort(candidate); % to find the closest, sort the candidates in order
    candidateLocation = candidateLocation(index); % save the location

    nLow = sum(repmat(candidate',length(vKprime),1)<vKprime,2); % using the sorted vector, find the period where the capital stock is closest to vKprime from below
    nLow(nLow<=1) = 1; % the location cannot go below 1.
    nLow(nLow>=length(index)) = length(index)-1; % the location cannot go over the length(index)-1: note that it's the closest from BELOW.
    nHigh = nLow+1; %define the period where the capital stock is closest to vKprime from above
    weightLow = (candidate(nHigh) - vKprime)./(candidate(nHigh)-candidate(nLow)); %compute the weight on the lower side
    weightLow(weightLow<0) = 0; % optional restriction on the extrapolation
    weightLow(weightLow>1) = 1; % optional restriction on the extrapolation
    
    LLow    = ((1-pAalpha)*Aprime*vKprime.^pAalpha./(pEeta*vC(candidateLocation(nLow )).^pRiskAversion)).^(pFrisch/(1+pAalpha*pFrisch));
    LHigh   = ((1-pAalpha)*Aprime*vKprime.^pAalpha./(pEeta*vC(candidateLocation(nHigh )).^pRiskAversion)).^(pFrisch/(1+pAalpha*pFrisch));
    rLow    = pAalpha*Aprime*vKprime.^(pAalpha-1).*LLow.^(1-pAalpha) - pDdelta;
    rHigh   = pAalpha*Aprime*vKprime.^(pAalpha-1).*LHigh.^(1-pAalpha) - pDdelta;

    tempV1 = tempV1 + (futureShock ~= iAprime).* pBbeta .*...
                mTransA(iA,iAprime).*...
                (weightLow.*((1./vC(candidateLocation(nLow ))).^(pRiskAversion).*(1+rLow ) - (1-pDdelta)*vLambda(candidateLocation(nLow )))...
           + (1-weightLow).*((1./vC(candidateLocation(nHigh))).^(pRiskAversion).*(1+rHigh) - (1-pDdelta)*vLambda(candidateLocation(nHigh))) );

end

% for the realized future shock level on the simulated path
Lfuture    = ((1-pAalpha).*vGridA(futureShock).*vKprime.^pAalpha./(pEeta*vC(iFuture).^pRiskAversion)).^(pFrisch/(1+pAalpha*pFrisch));
rfuture    = pAalpha*vGridA(futureShock).*vKprime.^(pAalpha-1).*Lfuture.^(1-pAalpha) - pDdelta;
tempV1     = tempV1 + pBbeta*...
                mTransRealized.*((1./vC(iFuture)).^(pRiskAversion).*(1+rfuture) - (1-pDdelta)*vLambda(iFuture)); 

% update the allocations
tempC      = (1./(tempV1 + vLambda)).^(1/pRiskAversion);
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