%% A canonical New Keynesian model (Rotemberg)
% 2023.09.25
% Hanbaek Lee (hanbaeklee1@gmail.com)
% When you use the code, please cite the paper 
% "A Dynamically Consistent Global Nonlinear Solution 
% Method in the Sequence Space and Applications."
%=========================    
% this file is to compute the dsge allocations over a simulated path.
%=========================    
%=========================
% housekeeping
%=========================
clc;
clear variables;
close all;
fnpath = '../functions';
addpath(fnpath);

%=========================
% load ss
%=========================

dir = '../solutions/nkrotem_ss.mat';
load(dir);
ss = load('../solutions/nkrotem_ss.mat');

%=========================
% exogenous shock setup
%=========================

% TFP shock

pNumGridA = 5;
pPersistence = 0.95;
pVol = 0.009; % 0.012;
[mTransA, vGridA] = ...
fnTauchen(pPersistence, 0, pVol^2, pNumGridA, 2);
vGridA = exp(vGridA);


% preference shock

pNumGridBS  = 5;
pPersistence = 0.90; 
pVol = 0.0180; % 0.0124; % 0.016; 
[mTransBS, vGridBS] = ...
fnTauchen(pPersistence, 0, pVol^2, pNumGridBS, 2);
vGridBS = exp(vGridBS);
% vGridBS = 1;


% MP shock

pNumGridMP  = 3;
pPersistence = 0.0;
pVol = 0.0025; % 0.0025;
[mTransMP, vGridMP] = ...
fnTauchen(pPersistence, 0, pVol^2, pNumGridMP, 1);

% if you wish to shut down the monetary policy shock:
% pNumGridMP  = 1;
% vGridMP     = 0;
% mTransMP    = 1;


% combine

vGridA      = kron(kron(vGridA,ones(pNumGridBS,1)),ones(pNumGridMP,1));
vGridBS     = kron(kron(ones(pNumGridA,1),vGridBS),ones(pNumGridMP,1));
vGridMP     = kron(kron(ones(pNumGridA,1),ones(pNumGridBS,1)),vGridMP);
mTransA     = kron(kron(mTransA,mTransBS),mTransMP);
pNumGridA   = pNumGridA*pNumGridBS*pNumGridMP;


%=========================
%simulation path
%=========================
seed = 100;
rng(seed);  
% T = 5001;
T = 10001;
% T = 20001;
% T = 30001;
% T = 50001;
% T = 3001;
BURNIN = 1000;
requiredTime = T+BURNIN;
pInitialPoint = 1;
vSimPath = fnSimulator(pInitialPoint,mTransA,BURNIN+T);
vLSimPath = [vSimPath(1);vSimPath(1:(end-1))];

%=========================    
%expected allocation path
%=========================
vPiAgg  = ss.piAgg*ones(requiredTime,1);
vmc     = ss.mc*ones(requiredTime,1);
vn      = ss.n*ones(requiredTime,1);
vw      = ss.w*ones(requiredTime,1);
vc      = ss.c*ones(requiredTime,1);
vY      = ss.Y*ones(requiredTime,1);
vi      = ss.i*ones(requiredTime,1)+ normrnd(0,0.0001,requiredTime,1);

vPiAggnew   = zeros(size(vPiAgg));
vmcnew      = zeros(size(vPiAgg));
vwnew       = zeros(size(vPiAgg));
vnnew       = zeros(size(vPiAgg));
vcnew       = zeros(size(vPiAgg));
vYnew       = zeros(size(vPiAgg));
vinew       = zeros(size(vPiAgg));

%=========================    
%resume from the last one
%=========================    

dir = ['../solutions/WIP_nkrotem_bc.mat'];
% load(dir);


%=========================    
%preparation
%=========================    

%for GE loop and acceleration

error2      = 10;
tol_ge      = 1e-8;

            
% updating weight

weightOld   = 0.9900;%0.9990;
weightOld2  = 0.9900;


tic;
%%
%=========================
%RTM
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

%iteration preparation    
pNumIter = 1;  

while error2>tol_ge
%=========================
%backward solution
%=========================
vSuffStat = [ss.i;vi(1:end-1)];
vSuffStatPrime = vi;

tempVf = 0;
tempVi = 0;
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

    tempVf = tempVf + (futureShock ~= iAprime).*pBeta.*(vGridBS(iAprime)./vBSshock).*pPsi.*...
                mTransA(iA,iAprime).* (...
                weightLow.* ...
                ((vc(candidateLocation(nLow))./vc)).^(-pRiskAversion) ...
                .* (1+vPiAgg(candidateLocation(nLow))).*(vPiAgg(candidateLocation(nLow))-pPiTarget) ...
                .* (vY(candidateLocation(nLow))./vY) ...
                +  ...
                (1-weightLow).* ...
                ((vc(candidateLocation(nHigh))./vc)).^(-pRiskAversion) ...
                .* (1+vPiAgg(candidateLocation(nHigh))).*(vPiAgg(candidateLocation(nHigh))-pPiTarget) ...
                .* (vY(candidateLocation(nHigh))./vY)...
                );

    tempVi = tempVi + (futureShock ~= iAprime).*pBeta.*(vGridBS(iAprime)./vBSshock).*...
                mTransA(iA,iAprime).*(weightLow.* ...
                (vc(candidateLocation(nLow))./vc).^(-pRiskAversion).*(1./(1+vPiAgg(candidateLocation(nLow)))) ...
                + (1-weightLow).* ...
                (vc(candidateLocation(nHigh))./vc).^(-pRiskAversion).*(1./(1+vPiAgg(candidateLocation(nHigh)))) ...
                );    

    tempVc = tempVc + (futureShock ~= iAprime).*pBeta.*(vGridBS(iAprime)./vBSshock).*...
                mTransA(iA,iAprime).*(weightLow.* ...
                (vc(candidateLocation(nLow))).^(-pRiskAversion).*((1+vi)./(1+vPiAgg(candidateLocation(nLow)))) ...
                + (1-weightLow).* ...
                (vc(candidateLocation(nHigh))).^(-pRiskAversion).*((1+vi)./(1+vPiAgg(candidateLocation(nHigh)))) ...
                );     

end

tempVf = tempVf + pBeta.*(vGridBS(futureShock)./vBSshock).*pPsi.*...
                mTransException.*((vc(iFuture)./vc)).^(-pRiskAversion) ...
                .* (1+vPiAgg(iFuture)).*(vPiAgg(iFuture)-pPiTarget) ...
                .* (vY(iFuture)./vY);
tempVi = tempVi + pBeta.*(vGridBS(futureShock)./vBSshock).*...
                mTransException.*(vc(iFuture)./vc).^(-pRiskAversion).*(1./(1+vPiAgg(iFuture)));
tempVc = tempVc + pBeta.*(vGridBS(futureShock)./vBSshock).*...
                mTransException.*(vc(iFuture)).^(-pRiskAversion).*((1+vi)./(1+vPiAgg(iFuture)));

vinew = (1./tempVi)-1;
vmcnew = ((pElasticity-1) + pPsi*(1+vPiAgg).*(vPiAgg-pPiTarget) - tempVf)./pElasticity;
vwnew = vmcnew.*vA;
vnnew = (vwnew./((pEta.*vc.^pRiskAversion))).^pFrisch;
vYnew = vA.*vnnew;
vcnew = vYnew.*(1-(pPsi/2).*(vPiAgg-pPiTarget).^2);


%=========================    
% simulate forward
%=========================    
vYf = (((1/pEta)*(pElasticity-1)/pElasticity)^(pFrisch/(1+pFrisch*pRiskAversion))) ...
    .* (vA).^((1+pFrisch)/(1+pFrisch*pRiskAversion));

vipast      = [ss.i;vinew(1:end-1)];
vPiAggnew = ( ((1+vi)./((1+vipast).^pRhor)).^(1/(1-pRhor)).*pBeta.*(1+pPiTarget).^pTaylorPi./((vY./vYf).^pTaylorY.*exp(vMPshock))...
            ).^(1/(1+pTaylorPi)) - 1;
% an alternative Taylor rule:
% vPiAggnew = ( (1+vi).*pBeta.*(1+pPiTarget).^pTaylorPi./((vY./vYf).^pTaylorY.*exp(vMPshock)) ).^(1/(1+pTaylorPi))- 1;

% error
error2 = mean(([...
    vPiAgg  - vPiAggnew;...
    vmc     - vmcnew;...
    vn      - vnnew;...
    vw      - vwnew;...
    vc      - vcnew;...
    vY      - vYnew;...
    vi      - vinew;...    
    ]).^2);

vPiAgg  = weightOld*vPiAgg    + (1-weightOld)*vPiAggnew;
vmc     = weightOld*vmc       + (1-weightOld)*vmcnew;
vn      = weightOld*vn        + (1-weightOld)*vnnew;
vw      = weightOld*vw        + (1-weightOld)*vwnew;
vc      = weightOld*vc        + (1-weightOld)*vcnew;
vY      = weightOld*vY        + (1-weightOld)*vYnew;
vi      = weightOld*vi        + (1-weightOld)*vinew;

if (floor((pNumIter-1)/50) == (pNumIter-1)/50)
%=========================  
% Report
%========================= 
[~,vYhpfilter] = hpfilter(log(vY),1600);
fprintf(' \n');
fprintf('Market Clearing Results \n');
fprintf('Error: %.13f \n', error2);
fprintf('Output volatility: %.13f \n', std(vYhpfilter)*100);
fprintf(' \n');

%=========================  
% Plot
%=========================  
subplot(2,4,1);
plot(1:requiredTime,vPiAgg(1:requiredTime));hold on;
plot(1:requiredTime,vPiAggnew(1:requiredTime),'-.');
xlim([1,requiredTime]);
hold off;
legend("Predicted \pi","Realized \pi","location","northeast");

subplot(2,4,2);
plot(1:requiredTime,vc(1:requiredTime));hold on;
plot(1:requiredTime,vcnew(1:requiredTime),'-.');
xlim([1,requiredTime]);
hold off;
legend("Predicted c","Realized c","location","northeast");

subplot(2,4,3);
plot(1:requiredTime,vw(1:requiredTime));hold on;
plot(1:requiredTime,vwnew(1:requiredTime),'-.');
xlim([1,requiredTime]);
hold off;
legend("Predicted w","Realized w","location","northeast");

subplot(2,4,4);
plot(1:requiredTime,vn(1:requiredTime));hold on;
plot(1:requiredTime,vnnew(1:requiredTime),'-.');
xlim([1,requiredTime]);
hold off;
legend("Predicted n","Realized n","location","northeast");

subplot(2,4,5);
plot(1:requiredTime,vmc(1:requiredTime));hold on;
plot(1:requiredTime,vmcnew(1:requiredTime),'-.');
xlim([1,requiredTime]);
hold off;
legend("Predicted mc","Realized mc","location","northeast");

subplot(2,4,6);
plot(1:requiredTime,vY(1:requiredTime));hold on;
plot(1:requiredTime,vYnew(1:requiredTime),'-.');
xlim([1,requiredTime]);
hold off;
legend("Predicted Y","Realized Y","location","northeast");

subplot(2,4,7);
plot(1:requiredTime,vi(1:requiredTime));hold on;
plot(1:requiredTime,vinew(1:requiredTime),'-.');
xlim([1,requiredTime]);
hold off;
legend("Predicted i","Realized i","location","northeast");

subplot(2,4,8);
scatter(1-vn,vPiAgg);
title("Global Phillips Curve");
pause(0.1);

%=========================  
% Save (mid)
%=========================  
dir = ['../solutions/WIP_nkrotem_bc.mat'];
save(dir);
toc;

end

pNumIter = pNumIter+1;
% error2 = 10e-20;

end

%=========================  
% Save (final)
%=========================  
dir = ['../solutions/nkrotem_bc.mat'];
save(dir);
