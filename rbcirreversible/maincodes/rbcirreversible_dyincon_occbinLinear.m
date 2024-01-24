%%%%
%=========================
%housekeeping
%=========================
% clc;
% clear variables;
% close all;
% fnpath = '../functions';
% addpath(fnpath);

%=========================
% load ss
%=========================
ss = load('../solutions/rbcirreversible_ss.mat');
load('../solutions/rbcirreversible.mat');
globalSol = load('../solutions/rbcirreversible.mat');
occbinSol = load('../solutions/rbcirreversible_occbin.mat');
gdsgeSol  = load('../solutions/rbcirreversible_gdsge.mat');

%%
%=========================
% shock path
%=========================
% define the shock path the consistent with the original simulation.

TFP = globalSol.vGridA(globalSol.vSimPath);
TFPprime = [TFP(2:end);TFP(end)];
erra = log(TFPprime) - globalSol.pPersistence*log(TFP);
erraGrid = unique(erra);
vErraPath = zeros(size(erra));
for iTrans = 1:length(vErraPath)
vErraPath(iTrans) = find(erraGrid==erra(iTrans));
end
mProb = zeros(length(erraGrid),length(erraGrid));
for iGrid = 1:length(erraGrid)
mProb(:,iGrid) = sum(erra==erraGrid(iGrid))/length(erra);
end

%%
%=========================
% load occbin
%=========================
%load the solution solved by gdsge algorithm.
%timing adjustment is necessary
load '../solutions/rbcirreversible_occbin.mat';
% vK = (1+k_p([1,1,1:requiredTime-2]))*ss.K;
% vC = (1+c_p([1,1,1:requiredTime-2]))*ss.C;
% vLambda = lambdak_p([1,1,1:requiredTime-2]);

% for linear solution
vK = (1+k_l([1,1,1:requiredTime-2]))*ss.K;
vC = (1+c_l([1,1,1:requiredTime-2]))*ss.C;
vLambda = lambdak_l([1,1,1:requiredTime-2]);

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

    nLow = sum(repmat(candidate',length(vKprime),1)<vKprime,2); % using the sorted vector, find the period where the capital stock is closest to vKprime from below
    nLow(nLow<=1) = 1; % the location cannot go below 1.
    nLow(nLow>=length(index)) = length(index)-1; % the location cannot go over the length(index)-1: note that it's the closest from BELOW.
    nHigh = nLow+1; %define the period where the capital stock is closest to vKprime from above
    weightLow = (candidate(nHigh) - vKprime)./(candidate(nHigh)-candidate(nLow)); %compute the weight on the lower side
    weightLow(weightLow<0) = 0; % optional restriction on the extrapolation
    weightLow(weightLow>1) = 1; % optional restriction on the extrapolation
    
    rprime  = pAlpha.*Aprime.*vKprime.^(pAlpha-1) - pDelta;
    Lambdaprime = weightLow.*vLambda(candidateLocation(nLow )) + (1-weightLow).*vLambda(candidateLocation(nHigh ));

    tempV1 = tempV1 + (futureShock ~= iAprime).* pBeta.*...
                mTransA(iA,iAprime).*...
                (weightLow.*(1./vC(candidateLocation(nLow ))).^(pRiskAversion).*(1+rprime ) ...
           + (1-weightLow).*(1./vC(candidateLocation(nHigh))).^(pRiskAversion).*(1+rprime) ...
           - (1-pDelta)*Lambdaprime);

end

% update the allocations
rfuture    = pAlpha*vGridA(futureShock).*vKprime.^(pAlpha-1) - pDelta;
tempV1     = tempV1 + pBeta*...
                mTransRealized.*((1./vC(iFuture)).^(pRiskAversion).*(1+rfuture) - (1-pDelta)*vLambda(iFuture)); 
vLambdanew = vC.^(-pRiskAversion) - tempV1;
tempC      = (1./(tempV1 + vLambda)).^(1/pRiskAversion);
vr         = pAlpha.*vA.*vK.^(pAlpha-1) - pDelta;
vY         = vA.*vK.^(pAlpha);
vI         = vY - tempC;
EE         = (tempC - vC)./vC;

% frictional setup; irreversibility
vI(vI<=pPhi*ss.I) = pPhi*ss.I;

%=========================    
% step 2: simulate forward
%=========================   
vKpast = [ss.K;vK(1:end-1)];
vIpast = [ss.K*pDelta;vI(1:end-1)];

vKnew = (1-pDelta)*vKpast + vIpast;
vCnew = vGridA(vSimPath).*vKnew.^(pAlpha) - vI;
vLambdanew(vI>pPhi*ss.I) = 0;

error2 = mean(([...
    vK      - vKnew;...
    vC      - vCnew;...
    vLambda - vLambdanew...
    ]).^2);

errorK = vK - vKnew;

%=========================  
% Report
%========================= 
fprintf(' \n');
fprintf(' \n');
fprintf('======================== \n');
fprintf('Linear solution report \n');
fprintf('======================== \n');
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

% save
save('../solutions/rbcirreversible_occbinlin_EE.mat','EE');

%%
%=========================  
% dynamic consistency report
%========================= 
fprintf('\n');
fprintf('======================== \n');
fprintf('dynamic consistency report for linear solution \n');
fprintf('======================== \n');
disp(['max absolute error (in pct. of steady-state K): ', num2str(100*max(abs(errorK))/ss.K),'%']);
disp(['root mean sqaured error (in pct. of steady-state K): ', num2str(100*sqrt(mean(errorK.^2))/ss.K),'%']);
disp(['max euler equation error: ', num2str(100*max(EE)),'%']);
disp(['root mean sqaured euler equation error: ', num2str(100*mean(EE.^2).^0.5),'%']);
fprintf('\n');

%%
%=========================  
% business cycle report
%========================= 
fprintf('======================== \n');
fprintf('Business cycle statistics for the raw time series (for report) \n');
fprintf('======================== \n');
fprintf('mean (investment): %.3f \n', mean(vI));
fprintf('mean (consumption): %.3f \n', mean(vC));
fprintf('------------------------ \n');
fprintf('st. dev. (investment): %.3f \n', std(vI));
fprintf('st. dev. (consumption): %.3f \n', std(vC));
fprintf('------------------------ \n');
fprintf('skewness (investment): %.3f \n', skewness(vI));
fprintf('skewness (consumption): %.3f \n', skewness(vC));
fprintf('------------------------ \n');
fprintf('Kurtosis (investment): %.3f \n', kurtosis(vI));
fprintf('Kurtosis (consumption): %.3f \n', kurtosis(vC));
fprintf('\n');

%%
%=========================  
% Distribution of dynamic inconsistency
%========================= 
figure;
hist(100*errorK/ss.K,100);
xlim([-1,1]);
xlabel("Dynamic consistency error (in pct. of steady-state K)")
location = ['../figures/err_hist_occbinlinear.pdf'];
saveas(gcf, location);

