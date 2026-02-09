%% A canonical DMP model with exogenous separation 
%April. 2023
%repeated transition method (RTM)
%=========================
%housekeeping
%=========================
clc;
clear variables;
close all;
fnpath = '../functions';
addpath(fnpath);

%=========================
%load ss
%=========================
dir = '../solutions/dmpexo_ss.mat';
load(dir);
ss = load('../solutions/dmpexo_ss.mat');

%=========================
%grid setup
%=========================
% KS(1998) type
% pNumGridA   = 2;
% mTransA     = [0.875,0.125;0.125,0.875];
% vGridA      = [0.99,1.01];

% Tauchen
pNumGridA = 11;
pPersistence = 0.983;
pVol = 0.0100;
[mTransA, vGridA] = ...
    fnTauchen(pPersistence, 0, pVol^2, pNumGridA, 3);
vGridA = exp(vGridA);

%=========================
%simulation path
%=========================
seed = 100;
rng(seed);    
T = 5001;
BURNIN = 100;
requiredTime = T+BURNIN;
pInitialPoint = 1;
vSimPath = fnSimulator(pInitialPoint,mTransA,BURNIN+T);
vLSimPath = [vSimPath(1);vSimPath(1:(end-1))];

%=========================    
%expected allocation path
%=========================
vn = ss.n*ones(requiredTime,1)+ normrnd(0,0.00001,requiredTime,1);
vq = ss.q*ones(requiredTime,1);
vv = ss.v*ones(requiredTime,1);
vu = ss.u*ones(requiredTime,1);
vc = ss.c*ones(requiredTime,1);
vw = ss.w*ones(requiredTime,1);

vnnew = zeros(size(vv));
vqnew = zeros(size(vv));
vvnew = zeros(size(vv));
vunew = zeros(size(vv));
vcnew = zeros(size(vv));
vwnew = zeros(size(vv));

%=========================    
%resume from the last one
%=========================    
% load '../solutions/WIPdmpexo_bc.mat';

%=========================    
%preparation
%=========================    
%for GE loop and acceleration
error2      = 10;
weightOld   = 0.9500;%0.985;
tol_ge      = 1e-9;

tic;
%%
%=========================
%RTM
%=========================
iA = vSimPath;
iFuture = [(2:requiredTime)';requiredTime];
futureShock = vSimPath(iFuture);
vA  = vGridA(iA);

mTransException = zeros(size(vn));
for iTrans = 1:length(vn)
mTransException(iTrans,1) = mTransA(iA(iTrans),futureShock(iTrans));
end

%iteration preparation    
pNumIter = 1;  

while error2>tol_ge
%=========================
%backward solution
%=========================
vSuffStat = [ss.n;vn(1:end-1)];
vSuffStatPrime = vn;
tempV1 = 0;
tempV2 = 0;

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
    
    tempV2 = tempV2 + (futureShock ~= iAprime).* pa.beta.*...
                mTransA(iA,iAprime).*((weightLow.*(vc./vc(candidateLocation(nLow))).^(pa.sigma) ...
                + (1-weightLow).*(vc./vc(candidateLocation(nHigh))).^(pa.sigma)));

end

tempV1 = tempV1 + pa.beta.*mTransException.*(vc./vc(iFuture)).^(pa.sigma)...
    .*(vGridA(futureShock)-vw(iFuture)+(1-pa.lambda)*pa.kappa./vq(iFuture));
tempV1 = max(1e-5,tempV1);
tempV2 = tempV2 + pa.beta.*mTransException.*(vc./vc(iFuture)).^(pa.sigma);
tempq = pa.kappa./tempV1;
vqnew = tempq;
vqnew(vqnew>1) =1;

vwnew = (1-pa.eta)*pa.b + pa.eta*(vA + pa.kappa.*(vv./vu));
vExpSDF = tempV2;

%=========================    
% simulate forward
%=========================    
vvnew   = (vqnew./pa.m).^(1/(-pa.xi1)).*(1-vn);

npast = [ss.n;vn(1:end-1)];
vpast = [ss.v;vvnew(1:end-1)];
qpast = [ss.q;vqnew(1:end-1)];
    
vnnew = (1-pa.lambda)*npast + vvnew.*vqnew;
vunew = 1-vnnew;
vcnew = vnnew.*vA - pa.kappa.*vvnew + (1-vnnew)*pa.b;

error2 = mean(abs([vvnew-vv;vunew-vu;vnnew-vn;vcnew-vc;vqnew-vq]));

vv = (1-weightOld)*vvnew+weightOld*vv;
vu = (1-weightOld)*vunew+weightOld*vu;
vn = (1-weightOld)*vnnew+weightOld*vn;
vc = (1-weightOld)*vcnew+weightOld*vc;
vq = (1-weightOld)*vqnew+weightOld*vq;
vw = (1-weightOld)*vwnew+weightOld*vw;

if (floor((pNumIter-1)/50) == (pNumIter-1)/50)

%=========================  
% Summary
%=========================  
u = vu;
%quarterly aggregation
numGroup = floor(length(u)/3);
group = kron(1:numGroup,[1,1,1]);
if mod(length(u),3)==1
group = [kron(1:numGroup,[1,1,1]),numGroup+1];
elseif mod(length(u),3)==2
group = [kron(1:numGroup,[1,1,1]),numGroup+1,numGroup+1];
end
uQ = accumarray(group',u)/3;
[~,ufiltered] = hpfilter(log(uQ),1600);
% disp(std(ufiltered)*100);

%=========================  
% Report
%=========================  
fprintf(' \n');
fprintf('Market Clearing Results \n');
fprintf('Error: %.10f \n', error2);
fprintf('Unemp. volatility: %.10f \n', std(ufiltered)*100);

%=========================  
% Plot
%=========================  
subplot(3,2,1);
plot(1:requiredTime,vv(1:requiredTime)./vu(1:requiredTime));hold on;
plot(1:requiredTime,vvnew(1:requiredTime)./vunew(1:requiredTime),'-.');
xlim([1,requiredTime]);
hold off;
legend("Predicted \theta","Realized \theta","location","northeast");

subplot(3,2,2);
plot(1:requiredTime,vn(1:requiredTime));hold on;
plot(1:requiredTime,vnnew(1:requiredTime),'-.');
xlim([1,requiredTime]);
hold off;
legend("Predicted n","Realized n","location","northeast");

subplot(3,2,3);
plot(1:requiredTime,vq(1:requiredTime));hold on;
plot(1:requiredTime,vqnew(1:requiredTime),'-.');
xlim([1,requiredTime]);
hold off;
legend("Predicted q","Realized q","location","northeast");

subplot(3,2,4);
plot(1:requiredTime,vv(1:requiredTime));hold on;
plot(1:requiredTime,vvnew(1:requiredTime),'-.');
xlim([1,requiredTime]);
hold off;
legend("Predicted v","Realized v","location","northeast");

subplot(3,2,5);
plot(1:requiredTime,vw(1:requiredTime));hold on;
plot(1:requiredTime,vwnew(1:requiredTime),'-.');
xlim([1,requiredTime]);
hold off;
legend("Predicted w","Realized w","location","northeast");

subplot(3,2,6);
plot(1:requiredTime,1-vn(1:requiredTime));hold on;
plot(1:requiredTime,1-vnnew(1:requiredTime),'-.');
xlim([1,requiredTime]);
hold off;
legend("Predicted unrate","Realized unrate","location","northeast");

pause(0.01);

% error2 = 1e-20;
%=========================  
% Save (mid)
%=========================  
save '../solutions/WIPdmpexo_bc.mat';
toc;

end

pNumIter = pNumIter+1;

end
%=========================  
% Save (final)
%=========================  
save '../solutions/dmpexo_bc.mat';
