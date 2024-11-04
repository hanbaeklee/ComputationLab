%% Solving the model in Aiyagari (1994)
% 2023.09.25
% Hanbaek Lee (hanbaeklee1@gmail.com)
%=========================    
% this file is to compute the stationary equilibrium.
%=========================    
%summary of methods used in the code:
%1.VFI: monotonicity + accelerator
%2.non-stochastic simulation (iteration or eigenvector)
%=========================
%housekeeping
%=========================
clc;
clear variables;
close all;
fnPath = '../functions';
addpath(fnPath);

%=========================
%parameter setup
%=========================
pAalpha     = 0.36; %from Aiyagari (1994)
pBbeta      = 0.96; %from Aiyagari (1994)
pDdelta     = 0.08; %from Aiyagari (1994)

%=========================
%grid setup
%=========================
%idiosyncratic income shock
pMu = 0.0;
pRho = 0.900;
pSigma = 0.0872;
pNumGridz   = 7;
[mTransz, vGridz] = ...
fnTauchen(pRho, pMu, pSigma^2, pNumGridz, 3);
vGridz = exp(vGridz');                

%wealth grid
pNumGrida   = 50;
pNumGrida2  = 100; %finer grid for interpolation

%finer grid near lower wealth (Maliar, Maliar, and Valli, 2010)
vGridamin   = 0;
vGridamax   = 150;
x           = linspace(0,0.5,pNumGrida);
x2          = linspace(0,0.5,pNumGrida2);
y           = x.^7/max(x.^7);
y2          = x2.^7/max(x2.^7);
vGrida      = vGridamin+(vGridamax-vGridamin)*y;
vGrida2     = vGridamin+(vGridamax-vGridamin)*y2;

%=========================
%stationary labor supply (exogenous)
%=========================
%specify the options (for lecture purpose)
optionStochastic1   = 1;
optionStochastic2   = 1;
verbose             = true;         % ge loop interim reports on/off
%=========================
%optionStochastic1: for stationary dist. of idiosyncratic shock
%optionStochastic2: for stationary dist. of eq. allocation
%=========================
%if==1: iteration method
%if==2: eigenvector method
%=========================
if optionStochastic1 == 1
%=========================
%option1: iteration method
vL = ones(length(mTransz(:,1)))/sum(ones(length(mTransz(:,1))));
error = 10;
while error>0.00000001
   vLprime =  mTransz'*vL;
   error = abs(vL-vLprime);
   vL = vLprime;
end
disp(vL);
%=========================
elseif optionStochastic1 == 2
%=========================
%option2: eigenvector method
[vL,~] = eigs(mTransz',1);
vL = vL/sum(vL);
disp(vL);
end
%=========================
supplyL = vGridz*vL;

%%
%=========================
%initial guess
%=========================
%initial guess
K           = 6;
%updating rule
weightOld   = 0.9900;
%initial guess of value function
mValue      = repmat(0.01*vGrida',1,pNumGridz);
%define the size of matrix
mValueNew   = zeros(size(mValue));
mPolc       = zeros(size(mValue));
mPolaprime  = zeros(size(mValue));
%define the size of the distribution
currentDist = ones(pNumGrida2,pNumGridz);
currentDist = currentDist/(pNumGrida2*pNumGridz); %normalize

%for GE loop and acceleration
error2      = 10;
tol_ge      = 1e-8;
AccInterval = 20;
AccStarting = 30;
pNumIter_GE = 1 ;

%loop for GE
tic;
while error2>tol_ge

%given K, all the prices are known.
r   = pAalpha*(K/supplyL)^(pAalpha-1)-pDdelta;
mmu = r+pDdelta;
w   = (1-pAalpha)*(K/supplyL)^(pAalpha);

%for VFI routine
error = 10;
pNumIter = 1;

%=========================
%optimization
%=========================
%loop for VFI
while error>1e-8
%parfor iz = 1:pNumGridz
for iz = 1:pNumGridz    

    z = vGridz(iz);

    %expected future value
    expVal = mValue*mTransz';
    expVal = expVal(:,iz);

    %to utilize the monotonicity
    minWealth = vGridamin;

    for ia = 1:pNumGrida

    %-------------non-acceleration-------------%
    if (floor((pNumIter-1)/AccInterval) == (pNumIter-1)/AccInterval || pNumIter<=AccStarting)

    a = vGrida(ia);
    budget = w*z + (1+r)*a;

    %optimal saving decision
    aprime = fnOptFAST(pBbeta,budget,vGrida,expVal,pNumGrida,minWealth);

    %borrowing constraint
    aprime(aprime<vGridamin) = vGridamin;

    %linear interpolation for off-the-grid aprime
    aLow = sum(vGrida<aprime);
    aLow(aLow<=1) = 1;
    aLow(aLow>=pNumGrida) = pNumGrida-1;
    aHigh = aLow+1;
    weightLow = (vGrida(aHigh) - aprime)/(vGrida(aHigh)-vGrida(aLow));
    value = weightLow*expVal(aLow)+(1-weightLow)*expVal(aHigh);

    c = budget - aprime;

    %update
    minWealth = aprime;
    mValueNew(ia,iz) = log(c)+pBbeta*value;
    mPolc(ia,iz) = c;
    mPolaprime(ia,iz) = aprime;

    %-------------acceleration-------------%
    else

    c = mPolc(ia,iz);
    aprime = mPolaprime(ia,iz);

    aLow = sum(vGrida<aprime);
    aLow(aLow<1) = 1;
    aLow(aLow>=pNumGrida) = pNumGrida-1;
    aHigh = aLow+1;
    weightLow = (vGrida(aHigh) - aprime)/(vGrida(aHigh)-vGrida(aLow));
    value = weightLow*expVal(aLow)+(1-weightLow)*expVal(aHigh);

    mValueNew(ia,iz) = log(c)+pBbeta*value;

    end
    end
end

%iteration routine
error = max(max(max(abs(mValueNew-mValue))));
mValue = mValueNew;

% if (floor((pNumIter-1)/30) == (pNumIter-1)/30)
% DISP=['Iteration is in process : ', num2str(pNumIter),' X ',num2str(error)];
% disp(DISP)
% toc;
% end

pNumIter = pNumIter+1;

end

%%
%=========================
%interpolation - wealth policy
%=========================
%define the size of interpolated policy function
mPolaprime2 = zeros(size(currentDist));
%interpolate
if pNumGrida == pNumGrida2
    mPolaprime2 = mPolaprime;
else
for iz = 1:pNumGridz
    mPolaprime2(:,iz) = interp1(vGrida,mPolaprime(:,iz),vGrida2,"linear","extrap");
end
end

%=========================
if optionStochastic2 == 1
%=========================
%option1: iteration method
%=========================
%simulation - stochastic
%=========================
error = 10;
while error>1e-8
%define the size of new distribution
nextDist = zeros(size(currentDist));

for iz = 1:pNumGridz

    for ia = 1:pNumGrida2

        a = vGrida2(ia);
        nexta = mPolaprime2(ia,iz);
        LB = sum(vGrida2<nexta);
        LB(LB<=0) = 1;
        LB(LB>=pNumGrida2) = pNumGrida2-1;
        UB = LB+1;
        weightLB = (vGrida2(UB) - nexta)/(vGrida2(UB)-vGrida2(LB));
        weightLB(weightLB<0) = 0;
        weightLB(weightLB>1) = 1;
        weightUB = 1-weightLB;

        mass = currentDist(ia,iz);

        for futureiz = 1:pNumGridz

            nextDist(LB,futureiz) = ...
                nextDist(LB,futureiz)...
                +mass*mTransz(iz,futureiz)*weightLB;

            nextDist(UB,futureiz) = ...
                nextDist(UB,futureiz)...
                +mass*mTransz(iz,futureiz)*weightUB;

        end

    end

end

%simulation routine;
error = max(abs(nextDist-currentDist),[],"all");
currentDist = nextDist;

end

%=========================
elseif optionStochastic2 == 2
%=========================
%option2: eigenvector method
%=========================
%simulation - nonstochastic
%=========================
%eigenvector method
mPolicy = zeros(pNumGrida2*pNumGridz,pNumGrida2*pNumGridz);

vLocationCombineda = kron(1:pNumGrida2,ones(1,pNumGridz));
vLocationCombinedz = kron(ones(size(vGrida2)),1:pNumGridz);

for iLocation = 1:pNumGridz*pNumGrida2

    ia = vLocationCombineda(iLocation);
    iz = vLocationCombinedz(iLocation);

    a = vGrida2(ia);
    nexta = mPolaprime2(ia,iz);
    LB = sum(vGrida2<nexta);
    LB(LB<=0) = 1;
    LB(LB>=pNumGrida2) = pNumGrida2-1;
    UB = LB+1;
    weightLB = (vGrida2(UB) - nexta)/(vGrida2(UB)-vGrida2(LB));
    weightLB(weightLB<0) = 0;
    weightLB(weightLB>1) = 1;
    weightUB = 1-weightLB;

    for izprime = 1:pNumGridz

        mPolicy(iLocation,:) = mPolicy(iLocation,:)+(vLocationCombineda==LB).*(vLocationCombinedz==izprime) * weightLB * mTransz(iz,izprime);
        mPolicy(iLocation,:) = mPolicy(iLocation,:)+(vLocationCombineda==UB).*(vLocationCombinedz==izprime) * weightUB * mTransz(iz,izprime);

    end

end

mPolicytrans = mPolicy';
% mPolicytrans = sparse(mPolicy');
[currentDist0,~] = eigs(mPolicytrans,1);%eig(mPolicy');
currentDist0 = currentDist0(:,1)/sum(currentDist0(:,1));
currentDist = zeros(pNumGrida2,pNumGridz);

for iLocation = 1:pNumGridz*pNumGrida2

    ia = vLocationCombineda(iLocation);
    iz = vLocationCombinedz(iLocation);

    currentDist(ia,iz) = currentDist0(vLocationCombineda==ia & vLocationCombinedz==iz);

end
currentDist(currentDist<0) = 0;

end

%endogenous aggregate allocation
marginalDista = sum(currentDist,2);
endoK = vGrida2*marginalDista;

%error and update
error2 = abs(endoK - K);
K = K.*weightOld+endoK.*(1-weightOld);

% report only spasmodically
if verbose == true && (floor((pNumIter_GE-1)/50) == (pNumIter_GE-1)/50) || error2<= tol_ge
%=========================  
% interim report
%=========================  

fprintf(' \n');
fprintf('market clearing results \n');
fprintf('max error: %.15f \n', error2);
fprintf('capital rent: %.15f \n', r);
fprintf('wage: %.15f \n', w);
fprintf('aggregate capital: %.15f \n', K);

% plot
close all;
figure;
plot(vGrida2,currentDist);
title("The wealth distributions for different labor endowments","fontsize",15)
saveas(gcf,'../figures/dist_ss.jpg');
hold off;

pause(0.01);
toc;

end

pNumIter_GE = pNumIter_GE + 1;

end
