%% An RBC model with asset price, irreversibility, and endogenous labor supply with infinite Frisch elasticity
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
ss = load('../solutions/rbcassetirrendolaborinftyfrisch_ss.mat');
globalSol = load('../solutions/rbcassetirrendolaborinftyfrisch_bc.mat');
load('../solutions/rbcassetirrendolaborinftyfrisch_bc.mat');

%%
%=========================
%backward solution
%=========================
% vectorize the shock-related paths
iA = tsimpath;
ifuture = [(2:pathlength)';pathlength];
futureshock = tsimpath(ifuture);
vA  = vgridA(iA);

% prior calculation of time-series of the transition probabilities to the realized
% aggregate shocks on the simulated path
mTransRealized = zeros(size(tk));
for iTrans = 1:length(tk)
mTransRealized(iTrans,1) = mtransA(iA(iTrans),futureshock(iTrans));
end

tempV1 = 0;
for iAprime = 1:pnumgridA

    Aprime = vgridA(iAprime);

    % find a period where the future shock realization is the same as
    % iAprime and the capital stock is closest to vKprime from below and above.
    candidate = tk(find(tsimpath==iAprime)); % iso-shock periods
    candidateLocation = find(tsimpath==iAprime); % iso-shock period locations
    candidate(candidateLocation>pathlength-burnin) = []; % last burnin periods cannot be a candidate
    candidate(candidateLocation<burnin) = [];  % initial burnin periods cannot be a candidate
    candidateLocation(candidateLocation>pathlength-burnin) = []; % last burnin periods cannot be a candidate
    candidateLocation(candidateLocation<burnin) = [];  % initial burnin periods cannot be a candidate
    [candidate,index] = sort(candidate); % to find the closest, sort the candidates in order
    candidateLocation = candidateLocation(index); % save the location

    klow = sum(repmat(candidate',length(tkprime),1)<tkprime,2); % using the sorted vector, find the period where the capital stock is closest to vKprime from below
    klow(klow<=1) = 1; % the location cannot go below 1.
    klow(klow>=length(index)) = length(index)-1; % the location cannot go over the length(index)-1: note that it's the closest from BELOW.
    khigh = klow+1; %define the period where the capital stock is closest to vKprime from above
    weightlow = (candidate(khigh) - tkprime)./(candidate(khigh)-candidate(klow)); %compute the weight on the lower side
    weightlow(weightlow<0) = 0; % optional restriction on the extrapolation
    weightlow(weightlow>1) = 1; % optional restriction on the extrapolation
   
    lambdaprime = weightlow.*tlambda(candidateLocation(klow )) + (1-weightlow).*tlambda(candidateLocation(khigh ));    
    wprime = weightlow.*tw(candidateLocation(klow )) + (1-weightlow).*tw(candidateLocation(khigh));
    rprime = weightlow.*tr(candidateLocation(klow )) + (1-weightlow).*tr(candidateLocation(khigh));

    tempV1 = tempV1 + (futureShock ~= iAprime).* pbeta.*...
                mtransA(iA,iAprime).*...
                (weightlow.*(1./tc(candidateLocation(klow ))).*(rprime) ...
           + (1-weightlow).*(1./tc(candidateLocation(khigh))).*(rprime) );

end

% for the realized future shock level on the simulated path
wfuture    = tw(ifuture);
rfuture    = tr(ifuture);

% update the allocations
tempV1     = tempV1 + pbeta*...
                mTransRealized.*(1./tc(ifuture)).*(rfuture); 
tempC      = ((1-tlambda)./tempV1);
EE         = (tempC - tc)./tc;

%%
ts1 = tk;%(tsimpath==3);
ts2 = EE;%(tsimpath==3);

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