%% A canonical RBC model
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
% load the stead-state equilibrium allocations
%=========================
load '../solutions/rbc_ss.mat';
ss = load('../solutions/rbc_ss.mat');

%=========================
% aggregate shock
%=========================
% Tauchen method
pNumGridA = 7;
pPersistence = 0.95;
pVol = 0.009;
[mTransA, vGridA] = ...
fnTauchen(pPersistence, 0, pVol^2, pNumGridA, 3);
vGridA = exp(vGridA);

%=========================
% simulation path
%=========================
seed = 100;
rng(seed);    
% T = 1001;
T = 2001;
% T = 3001;
% T = 5001;
% T = 10001;
BURNIN = 500; %initial burn-in period
requiredTime = T+BURNIN; %number of total periods
pInitialPoint = 1; 
vSimPath = fnSimulator(pInitialPoint,mTransA,BURNIN+T);

%=========================    
% initial guess for the allocation path
%=========================
vr      = ss.r*ones(requiredTime,1);
vw      = ss.w*ones(requiredTime,1);
vL      = ss.L*ones(requiredTime,1);
vC      = ss.C*ones(requiredTime,1);
vK      = ss.K*ones(requiredTime,1) + normrnd(0,0.0000001,requiredTime,1);
vY      = ss.Y*ones(requiredTime,1);
vT      = ss.T*ones(requiredTime,1);
vI      = ss.K*pDdelta*ones(requiredTime,1);

% separate paths to be updated iteratively
vCnew   = zeros(requiredTime,1);
vKnew   = ss.K*ones(requiredTime,1);

%=========================    
% resume from the last one
%=========================    
% use the following line if you wish to start from where you stopped
% before.
% load load './solutions/WIP_rbc_bc.mat';

%=========================    
% preparation
%=========================    
% the updating weights
weightOld1   = 0.9000; % updating weight for capital stock 
weightOld2   = 0.9000; % updating weight for consumption
weightOld3   = 0.9000; % updating weight for wage

% for an extremely high accuracy, you might consider weight as high as
% 0.9990;

tic;
%%
%=========================
% repeated transition method
%=========================
% iteration preparation    
pNumIter    = 1;  % this is for interim reports
error2      = 10; % this is for terminal condition

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

while error2>1e-8
    
%=========================
% step 1: backward solution
%=========================
% even if it's the bacward-solution step, it does not include the backward
% loop, as the step is vectorized.

% calculate the future endogenous capital allocations based on the time
% series of endogenous capital in the (n-1)th iteration.
vKprime = [vK(2:end);vK(1)];

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

    KLow = sum(repmat(candidate',length(vKprime),1)<vKprime,2); % using the sorted vector, find the period where the capital stock is closest to vKprime from below
    KLow(KLow<=1) = 1; % the location cannot go below 1.
    KLow(KLow>=length(index)) = length(index)-1; % the location cannot go over the length(index)-1: note that it's the closest from BELOW.
    KHigh = KLow+1; %define the period where the capital stock is closest to vKprime from above
    weightLow = (candidate(KHigh) - vKprime)./(candidate(KHigh)-candidate(KLow)); %compute the weight on the lower side
    weightLow(weightLow<0) = 0; % optional restriction on the extrapolation
    weightLow(weightLow>1) = 1; % optional restriction on the extrapolation
    
    LLow    = 1-(pTtheta.*vC(candidateLocation(KLow )))./(vw(candidateLocation(KLow )));
    LHigh   = 1-(pTtheta.*vC(candidateLocation(KHigh )))./(vw(candidateLocation(KHigh )));
    rLow    = pAalpha*Aprime*vKprime.^(pAalpha-1).*LLow.^(1-pAalpha) - pDdelta;
    rHigh   = pAalpha*Aprime*vKprime.^(pAalpha-1).*LHigh.^(1-pAalpha) - pDdelta;

    tempV1 = tempV1 + (futureShock ~= iAprime).* pBbeta.*...
                mTransA(iA,iAprime).*...
                (weightLow.*(1./vC(candidateLocation(KLow ))).^(pRiskAversion).*(1+rLow ) ...
           + (1-weightLow).*(1./vC(candidateLocation(KHigh))).^(pRiskAversion).*(1+rHigh) ...
           );

end

% for the realized future shock level on the simulated path
Lfuture     = 1-(pTtheta.*vC(iFuture))./(vw(iFuture));
rfuture     = pAalpha*vGridA(futureShock).*vKprime.^(pAalpha-1).*Lfuture.^(1-pAalpha) - pDdelta;
tempV1      = tempV1 + pBbeta*...
              mTransRealized.*((1./vC(iFuture)).^(pRiskAversion).*(1+rfuture)); 

% update the allocations
tempC       = ((1./tempV1).^(1/pRiskAversion));
vL          = 1-(pTtheta.*tempC)./vw;
vr          = pAalpha*vA.*vK.^(pAalpha-1).*vL.^(1-pAalpha) - pDdelta;
vY          = vA.*vK.^(pAalpha).*vL.^(1-pAalpha);
vI          = vY - tempC;
vwnew       = (1-pAalpha)*vA.*vK.^(pAalpha).*vL.^(-pAalpha);

%=========================    
% step 2: simulate forward
%=========================   
vKpast = [vK(end);vK(1:end-1)];
vIpast = [vI(end);vI(1:end-1)];

vKnew = (1-pDdelta)*vKpast + vIpast;
vCnew = vGridA(vSimPath).*vKnew.^(pAalpha).*vL.^(1-pAalpha) - vI;

error2 = mean(([...
    vC      - vCnew;...
    vK      - vKnew...
    ]).^2);
errorK = vK - vKnew;

vC      = weightOld1*vC        + (1-weightOld1)*vCnew;
vK      = weightOld2*vK        + (1-weightOld2)*vKnew;
vw      = weightOld3*vw        + (1-weightOld3)*vwnew;

if (floor((pNumIter-1)/50) == (pNumIter-1)/50)
%=========================  
% Report
%========================= 
Phrase = ['Iteration is in progress: ',num2str(pNumIter),'st iteration'];
disp(Phrase);
fprintf('Convergence criterion: \n');
fprintf('Error: %.18f \n', error2);
fprintf(' \n');
    
subplot(1,2,1);
plot(1:requiredTime,vK(1:requiredTime));hold on;
plot(1:requiredTime,vKnew(1:requiredTime),'-.');
xlim([1,requiredTime]);
hold off;
legend("Predicted K","Realized K","location","northeast");

subplot(1,2,2);
plot(1:requiredTime,vC(1:requiredTime));hold on;
plot(1:requiredTime,vCnew(1:requiredTime),'-.');
xlim([1,requiredTime]);
hold off;
legend("Predicted C","Realized C","location","northeast");

pause(0.2);

%=========================  
% save (mid)
%=========================  
save '../solutions/WIP_rbc_bc.mat';
% toc;

end

pNumIter = pNumIter+1;

end % end of the final loop
toc;

%=========================  
% save (final)
%=========================  
save '../solutions/rbc_bc.mat';

%%
%=========================  
% final report
%========================= 
fprintf('\n');
fprintf('======================== \n');
fprintf('Final report\n');
fprintf('======================== \n');
fprintf('Convergence criterion: \n');
fprintf('Error: %.9f \n', error2);

fprintf('\n');
fprintf('======================== \n');
fprintf('Business cycle statistics for the raw time series\n');
fprintf('======================== \n');
fprintf('mean log(output): %.4f \n', mean(log(vY)));
fprintf('st. dev. log(output): %.4f \n', std(log(vY)));
fprintf('skewness log(output): %.4f \n', skewness(log(vY)));
fprintf('------------------------ \n');
fprintf('mean log(investment): %.4f \n', mean(log(vI)));
fprintf('st. dev. log(investment): %.4f \n', std(log(vI)));
fprintf('skewness log(investment): %.4f \n', skewness(log(vI)));
fprintf('------------------------ \n');
fprintf('mean log(consumption): %.4f \n', mean(log(vC)));
fprintf('st. dev. log(consumption): %.4f \n', std(log(vC)));
fprintf('skewness log(consumption): %.4f \n', skewness(log(vC)));

fprintf('\n');
fprintf('======================== \n');
fprintf('Business cycle statistics for the HP-filtered time series\n');
fprintf('======================== \n');
[~,vYhpfilter] = hpfilter(log(vY),1600);
[~,vIhpfilter] = hpfilter(log(vI),1600);
[~,vChpfilter] = hpfilter(log(vC),1600);
fprintf('mean log(output): %.4f \n', mean(log(vYhpfilter)));
fprintf('st. dev. log(output): %.4f \n', std(log(vYhpfilter)));
fprintf('skewness log(output): %.4f \n', skewness(log(vYhpfilter)));
fprintf('------------------------ \n');
fprintf('mean log(investment): %.4f \n', mean(log(vIhpfilter)));
fprintf('st. dev. log(investment): %.4f \n', std(log(vIhpfilter)));
fprintf('skewness log(investment): %.4f \n', skewness(log(vIhpfilter)));
fprintf('------------------------ \n');
fprintf('mean log(consumption): %.4f \n', mean(log(vChpfilter)));
fprintf('st. dev. log(consumption): %.4f \n', std(log(vChpfilter)));
fprintf('skewness log(consumption): %.4f \n', skewness(log(vChpfilter)));

%%
%=========================  
% dynamic consistency report
%========================= 
fprintf('\n');
fprintf('======================== \n');
fprintf('dynamic consistency report for rtm \n');
fprintf('======================== \n');
disp(['max absolute error (in pct. of steady-state K): ', num2str(100*max(abs(errorK))/ss.K),'%']);
disp(['root mean sqaured error (in pct. of steady-state K): ', num2str(100*sqrt(mean(errorK.^2))/ss.K),'%']);
fprintf('\n');
figure;
hist(100*errorK/ss.K,100);
xlim([-1,1]);
xlabel("Dynamic consistency error (in pct. of steady-state K)")
location = ['../figures/err_hist.pdf'];
saveas(gcf, location);

%%
%=========================  
% fitting LoM into the linear specification 
%========================= 
endoState = vK(BURNIN+1:end-1);
exoState = vGridA(vSimPath(BURNIN+1:end-1));
endoStatePrime = vK(BURNIN+2:end);
intercept = ones(size(endoState));

% independent variable
x = [intercept ...
    ,log(endoState) ...
    ,log(exoState) ...
    ,log(endoState).*log(exoState)...
    ];

% dependent variable
y = log(endoStatePrime);

[coeff,bint1,r1,rint1,R1] = regress(y,x);
fprintf('======================== \n');
fprintf('Fitting the true LoM into the log-linear specification\n');
fprintf('======================== \n');
disp(['R-squared: ',num2str(R1(1))]);
fprintf(' \n');
fitlm(x,y,'Intercept',false)

% recover the implied dynamics
startingPoint = BURNIN+1;
startingEndo = endoState(1);
recovered = ones(1,requiredTime - BURNIN)*startingEndo;
for iTrans = 1:(requiredTime - BURNIN-1)
% endoStateTemp = endoState(iTrans);
endoStateTemp = recovered(iTrans);
exoStateTemp = exoState(iTrans);
tempX = [1 ...
    ,log(endoStateTemp)...
    ,log(exoStateTemp)...
    ,log(endoStateTemp)*log(exoStateTemp)...
    ];
tempVal = coeff'*tempX';
recovered(iTrans+1) = exp(tempVal);
end

samplePeriod = 500:1000;
figure;
plot(samplePeriod,endoState(samplePeriod),'Color','red','LineWidth',1.5);hold on;
plot(samplePeriod,recovered(samplePeriod),'Color','blue','LineStyle','--','LineWidth',1.5);
legend("True LoM","Linear LoM","location","best","FontSize",15)
location = ['../figures/lom.pdf'];
saveas(gcf, location);


%%
%=========================  
% fitting wage dynamics into the linear specification 
%========================= 
endostate = vK(BURNIN+1:end-1);
exostate = vGridA(vSimPath(BURNIN+1:end-1));
endoprice = vw(BURNIN+1:end-1);
intercept = ones(size(endostate));

% independent variable
x = [intercept ...
    ,log(endostate) ...
    ,log(exostate) ...
    ,log(endostate).*log(exostate)...
    ];

% dependent variable
y = log(endoprice);

[coeff,bint1,r1,rint1,R1] = regress(y,x);
fprintf('======================== \n');
fprintf('Fitting the true LoM into the log-linear specification\n');
fprintf('======================== \n');
disp(['R-squared: ',num2str(R1(1))]);
fprintf(' \n');
fitlm(x,y,'Intercept',false)

% recover the implied dynamics
recovered = zeros(1,requiredTime - BURNIN);
for itrans = 1:(requiredTime - BURNIN-1)
% endoStateTemp = endoState(iTrans);
endostatetemp = endostate(itrans);
exostatetemp = exostate(itrans);
tempX = [1 ...
    ,log(endostatetemp)...
    ,log(exostatetemp)...
    ,log(endostatetemp)*log(exostatetemp)...
    ];
tempVal = coeff'*tempX';
recovered(itrans) = exp(tempVal);
end

samplePeriod = 500:1000;
figure;
plot(samplePeriod,endoprice(samplePeriod),'Color','red','LineWidth',1.5);hold on;
plot(samplePeriod,recovered(samplePeriod),'Color','blue','LineStyle','--','LineWidth',1.5);
legend("True LoM","Linear LoM","location","best","FontSize",15)
location = ['../figures/lom_w.pdf'];
saveas(gcf, location);

