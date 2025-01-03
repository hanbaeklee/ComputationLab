%% A canonical New Keynesian model (Rotemberg)
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
clc;
clear variables;
close all; 
fnPath = './functions';
addpath(fnPath);

%=========================
%load ss
%=========================
ss = load('../solutions/nkrotem_ss.mat');
globalSol = load('../solutions/nkrotem_bc.mat');
load('../solutions/nkrotem_bc.mat');

%%
%=========================
%backward solution
%=========================
iA = vSimPath;
iFuture = [(2:requiredTime)';requiredTime];
futureShock = vSimPath(iFuture);
vA  = vGridA(iA);
vBSshock = vGridBS(iA);
vMPshock = vGridMP(iA);

mTransException = zeros(size(vPiAgg));
for iTrans = 1:length(vPiAgg)
mTransException(iTrans,1) = mTransA(iA(iTrans),futureShock(iTrans));
end

vSuffStat = [ss.i;vi(1:end-1)];
vSuffStatPrime = vi;

tempVc = 0;
for iAprime = 1:pNumGridA

    Aprime = vGridA(iAprime);

    % find a period where the future shock realization is the same as
    % iAprime and the sufficient stat that is closest to target stat from below and above.
    candidate = vSuffStat(find(vSimPath==iAprime)); % iso-shock periods
    candidateLocation = find(vSimPath==iAprime); % iso-shock period locations
    candidate(candidateLocation>requiredTime-BURNIN) = []; % last burnin periods cannot be a candidate
    candidate(candidateLocation<BURNIN) = []; % initial burnin periods cannot be a candidate
    candidateLocation(candidateLocation>requiredTime-BURNIN) = []; % last burnin periods cannot be a candidate
    candidateLocation(candidateLocation<BURNIN) = []; % initial burnin periods cannot be a candidate
    [candidate,index] = sort(candidate); % to find the closest, sort the candidates in order
    candidateLocation = candidateLocation(index); % save the location

    nLow = sum(repmat(candidate',length(vSuffStatPrime),1)<vSuffStatPrime,2); % using the sorted vector, find the period where the sufficient stat is closest to the target from below
    nLow(nLow<=1) = 1; % the location cannot go below 1.
    nLow(nLow>=length(index)) = length(index)-1; % the location cannot go over the length(index)-1: note that it's the closest from BELOW.
    nHigh = nLow+1; % define the period where the sufficient stat is closest to the taget stat from above
    weightLow = (candidate(nHigh) - vSuffStatPrime)./(candidate(nHigh)-candidate(nLow)); %compute the weight on the lower side
    weightLow(weightLow<0) = 0; % optional restriction on the extrapolation
    weightLow(weightLow>1) = 1; % optional restriction on the extrapolation
    
    weightLow(candidate(nHigh)==candidate(nLow)) = 1;

    tempVc = tempVc + (futureShock ~= iAprime).*pBeta.*(vGridBS(iAprime)./vBSshock).*...
            mTransA(iA,iAprime).*(weightLow.* ...
            (vc(candidateLocation(nLow))).^(-pRiskAversion).*((1+vi)./(1+vPiAgg(candidateLocation(nLow)))) ...
            + (1-weightLow).* ...
            (vc(candidateLocation(nHigh))).^(-pRiskAversion).*((1+vi)./(1+vPiAgg(candidateLocation(nHigh)))) ...
            );   

end

% for the realized future shock level on the simulated path
tempVc = tempVc + pBeta.*(vGridBS(futureShock)./vBSshock).*...
                mTransException.*(vc(iFuture)).^(-pRiskAversion).*((1+vi)./(1+vPiAgg(iFuture)));


% update the allocations
tempC       = ((1./tempVc).^(1/pRiskAversion));
EE         = (tempC - vc)./vc;

%%
ts1 = vSuffStat;%(vSimPath==3);
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
xlabel("Nominal interest rate","FontSize",15);
ylabel("Euler equation error (in Log)","FontSize",15);
location = ['../figures/ee.pdf'];
saveas(gcf, location);