%% A canonical DMP model with exogenous separation 
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
ss = load('../solutions/dmpexo_ss.mat');
globalSol = load('../solutions/dmpexo_bc.mat');
load('../solutions/dmpexo_bc.mat');

%%
%=========================
%backward solution
%=========================
iA = vSimPath;
iFuture = [(2:requiredTime)';requiredTime];
futureShock = vSimPath(iFuture);
vA  = vGridA(iA);

mTransException = zeros(size(vn));
for iTrans = 1:length(vn)
mTransException(iTrans,1) = mTransA(iA(iTrans),futureShock(iTrans));
end

vSuffStat = [ss.n;vn(1:end-1)];
vSuffStatPrime = vn;
tempV1 = 0;
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
       
    tempV1 = tempV1 + (futureShock ~= iAprime).* pa.beta.*...
                mTransA(iA,iAprime).*(...
                  (weightLow.*(vc./vc(candidateLocation(nLow))).^(pa.sigma)...
                .* (Aprime-vw(candidateLocation(nLow)) + (1-pa.lambda)*pa.kappa./vq(candidateLocation(nLow))) ...
                + (1-weightLow).*(vc./vc(candidateLocation(nHigh))).^(pa.sigma)...
                .* (Aprime-vw(candidateLocation(nHigh)) + (1-pa.lambda)*pa.kappa./vq(candidateLocation(nHigh)))));
    
end

tempV1 = tempV1 + pa.beta.*mTransException.*(vc./vc(iFuture)).^(pa.sigma)...
    .*(vGridA(futureShock)-vw(iFuture)+(1-pa.lambda)*pa.kappa./vq(iFuture));
tempq = pa.kappa./tempV1;
tempq(tempq>1) =1;

% update the allocations
EE         = (tempq - vq)./vq;

%%
ts1 = vSuffStat(BURNIN+1:end-BURNIN);%(vSimPath==3);
ts2 = EE(BURNIN+1:end-BURNIN);%(vSimPath==3);

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
xlabel("Sufficient statistics","FontSize",15);
ylabel("Euler equation error (in Log)","FontSize",15);
location = ['../figures/ee.pdf'];
saveas(gcf, location);