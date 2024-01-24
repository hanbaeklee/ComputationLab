%% RBC model with irreversible investment
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
ss = load('../solutions/rbcassetprice_ss.mat');
globalSol = load('../solutions/rbcassetprice_bc.mat');
load('../solutions/rbcassetprice_bc.mat');

%%
%=========================
%backward solution
%=========================
% vectorize the shock-related paths
iA = vSimPath;
iFuture = [(2:requiredTime)';requiredTime];
futureShock = vSimPath(iFuture);
vA  = vGridA(iA);

% prior calculation of time-series of the transition probabilities to the realized
% aggregate shocks on the simulated path
mTransRealized = zeros(size(vK));
for iTrans = 1:length(vK)
mTransRealized(iTrans,1) = mTransA(iA(iTrans),futureShock(iTrans));
end

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

    KLow = sum(repmat(candidate',length(vKprime),1)<vKprime,2); % using the sorted vector, find the period where the capital stock is closest to vKprime from below
    KLow(KLow<=1) = 1; % the location cannot go below 1.
    KLow(KLow>=length(index)) = length(index)-1; % the location cannot go over the length(index)-1: note that it's the closest from BELOW.
    KHigh = KLow+1; %define the period where the capital stock is closest to vKprime from above
    weightLow = (candidate(KHigh) - vKprime)./(candidate(KHigh)-candidate(KLow)); %compute the weight on the lower side
    weightLow(weightLow<0) = 0; % optional restriction on the extrapolation
    weightLow(weightLow>1) = 1; % optional restriction on the extrapolation
   
    rprime  = pAalpha*Aprime*vKprime.^(pAalpha-1) + (1-pDdelta);

    tempV1 = tempV1 + (futureShock ~= iAprime).* pBbeta.*...
                mTransA(iA,iAprime).*...
                (weightLow.*(1./vC(candidateLocation(KLow ))).^(pRiskAversion).*(rprime) ...
           + (1-weightLow).*(1./vC(candidateLocation(KHigh))).^(pRiskAversion).*(rprime) );

end

% for the realized future shock level on the simulated path
rfuture    = pAalpha.*vGridA(futureShock).*vKprime.^(pAalpha-1) + (1-pDdelta);
tempV1     = tempV1 + pBbeta*...
                mTransRealized.*(1./vC(iFuture)).^(pRiskAversion).*(rfuture); 

% update the allocations
tempC      = (1./tempV1).^(1/pRiskAversion);
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