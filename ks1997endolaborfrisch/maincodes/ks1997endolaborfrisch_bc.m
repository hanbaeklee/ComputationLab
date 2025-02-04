%% Krusell and Smith (1997) with endogenous labor supply
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
% macro choices
%=========================
verbose1        = false;        % vfi interim reports on/off
verbose2        = true;         % ge loop interim reports on/off

%=========================
% load the stationary equilibrium allocations
%=========================
load '../solutions/ks1997endolabor_ss.mat';
ss = load('../solutions/ks1997endolabor_ss.mat');
pBbar = 20;
%=========================
% numerical parameters
%=========================
tol_ge          = 1e-8;
weightold1      = 0.9500;
weightold2      = 0.9500;
weightold3      = 0.9500;
weightold4      = 0.9500;
weightold5      = 0.9990;

%=========================
% aggregate shock
%=========================
% aggregate productivity shock
pnumgridA   = 2;
mtransA     = [0.875,0.125;0.125,0.875];
vgridA      = [0.99,1.01];
% vgridA      = [1.00,1.00];

%=========================
% simulation path
%=========================
seed = 100;
rng(seed);
T = 1001;% T = 1001; % T = 3001;% T = 5001;% T = 10001;
burnin = 500;
pathlength = T+burnin;
pinitialpoint = 1;
tsimpath = fnSimulator(pinitialpoint,mtransA,burnin+T);

%=========================        
% declare equilibrium objects
%=========================    
mpolc       = zeros(pnumgridomega,pnumgridz,pathlength);
mpolkprime  = zeros(pnumgridomega,pnumgridz,pathlength);
mpolkprime_new = zeros(pnumgridomega,pnumgridz,pathlength);
mpolbprime  = zeros(pnumgridomega,pnumgridz,pathlength);
mpolbprime_new = zeros(pnumgridomega,pnumgridz,pathlength);
mlambda     = zeros(pnumgridomega,pnumgridz,pathlength);
mlambda_new = zeros(pnumgridomega,pnumgridz,pathlength);
mpolc_new   = zeros(pnumgridomega,pnumgridz,pathlength);
mpoln       = zeros(pnumgridomega,pnumgridz,pathlength);
mphi        = zeros(pnumgridomega,pnumgridz,pathlength);
mphi_new    = zeros(pnumgridomega,pnumgridz,pathlength);

%=========================            
% start and end points
%=========================    
startingdist = currentdist;
for itrans = 1:pathlength
mpolkprime(:,:,itrans) = ss.mpolkprime;
mpolbprime(:,:,itrans) = ss.mpolbprime;
mpoln(:,:,itrans) = ss.mpoln;
mlambda(:,:,itrans) = ss.mlambda;
mphi(:,:,itrans) = ss.mphi;
mpolc(:,:,itrans) = ss.mpolc;
end

%=========================     
% initial guess
%========================= 
tOmega = ss.Omega*ones(pathlength,1);
tK = ss.K*ones(pathlength+1,1)+ normrnd(0,0.0001,pathlength+1,1);
tB = ss.endoB*ones(pathlength+1,1);

% for the initial iteration, slightly perturbed capital path is necessary.
tsupplyL = ss.supplyL*ones(pathlength,1);
tY = zeros(pathlength,1);
tC = zeros(pathlength,1);
tq = ones(pathlength,1)*ss.q;
tlambda =  sum(ss.mlambda.*ss.currentdist,"all")*ones(pathlength,1);
tphi    =  sum(ss.mphi.*ss.currentdist,"all")*ones(pathlength,1);

%%
%=========================    
% resume from the last one
%=========================    
% use the following line if you wish to start from where you stopped
% before.
load '../solutions/WIP_ks1997endolabor_bc.mat';

%%        
%=========================    
% outer loop
%=========================    
tic;
error2 = 10;
pnumiter_ge = 1;

while error2>tol_ge
     
%=========================    
% 2. backward solution
%=========================    
for itrans = pathlength:-1:1

iA = tsimpath(itrans);
A = vgridA(iA);
K = tK(itrans);
supplyL = tsupplyL(itrans);

%given K, all the prices are known.
r   = palpha*A*(K/supplyL)^(palpha-1)-pdelta;
tr(itrans)  = r;
mmu = r+pdelta;
w   = (1-palpha)*A*(K/supplyL)^(palpha);
q   = tq(itrans);

if itrans == pathlength
futureShock = tsimpath(itrans);
ifuture = itrans;
Kprime = tK(itrans+1);
else
futureShock = tsimpath(itrans+1);
ifuture = itrans+1;
Kprime = tK(itrans+1);
end

% expected future value (rational expectation)
mexpectation = 0;
mexpectation2 = 0;
for iAprime = 1:pnumgridA
    Aprime = vgridA(iAprime);

    for izprime = 1:pnumgridz
    zprime = vgridz(izprime);
    
    if futureShock ~= iAprime || itrans == pathlength      
    % find a period where the future shock realization is the same as
    % iAprime and the capital stock is closest to Kprime from below and above.
    candidate = tK(find(tsimpath==iAprime)); % iso-shock periods
    candidatelocation = find(tsimpath==iAprime); % iso-shock period locations
    candidate(candidatelocation>pathlength-burnin) = []; % last burnin periods cannot be a candidate
    candidate(candidatelocation<burnin) = [];  % initial burnin periods cannot be a candidate
    candidatelocation(candidatelocation>pathlength-burnin) = []; % last burnin periods cannot be a candidate
    candidatelocation(candidatelocation<burnin) = [];  % initial burnin periods cannot be a candidate
    [candidate,index] = sort(candidate); % to find the closest, sort the candidates in order
    candidatelocation = candidatelocation(index); % save the location

    Klow = sum(candidate<Kprime); % using the sorted vector, find the period where the capital stock is closest to Kprime from below
    Klow(Klow<=1) = 1; % the location cannot go below 1.
    Klow(Klow>=length(index)) = length(index)-1; % the location cannot go over the length(index)-1: note that it's the closest from BELOW.
    Khigh = Klow+1; %define the period where the capital stock is closest to Kprime from above
    weightlow = (candidate(Khigh) - Kprime)./(candidate(Khigh)-candidate(Klow)); %compute the weight on the lower side
    weightlow(weightlow<0) = 0; % optional restriction on the extrapolation
    weightlow(weightlow>1) = 1; % optional restriction on the extrapolation
    
    K2Lprimelow = Kprime/tsupplyL(candidatelocation(Klow));
    rprimelow   = palpha*Aprime*(K2Lprimelow)^(palpha-1)-pdelta;
    wprimelow   = (1-palpha)*Aprime*(K2Lprimelow)^(palpha);        

    K2Lprimehigh = Kprime/tsupplyL(candidatelocation(Khigh));
    rprimehigh   = palpha*Aprime*(K2Lprimehigh)^(palpha-1)-pdelta;
    wprimehigh   = (1-palpha)*Aprime*(K2Lprimehigh)^(palpha);        

    rprime      = weightlow*rprimelow + (1-weightlow)*rprimehigh;
    wprime      = weightlow*wprimelow + (1-weightlow)*wprimehigh;
    qprime      = weightlow*tq(candidatelocation(Klow)) + (1-weightlow)*tq(candidatelocation(Khigh));

    mpolomegaprime = mpolkprime(:,:,itrans)*(1+rprime) + mpolbprime(:,:,itrans);

    mpolkprimetemp = weightlow*mpolkprime(:,izprime,candidatelocation(Klow)) ...
                + (1-weightlow)*mpolkprime(:,izprime,candidatelocation(Khigh));
    mpolkprimeprime = interp1(vgridomega',mpolkprimetemp,...
        mpolomegaprime,"linear","extrap");
    mpolkprimeprime(mpolkprimeprime<0) = 0;
    
    mpolbprimetemp = weightlow*mpolbprime(:,izprime,candidatelocation(Klow)) ...
                + (1-weightlow)*mpolbprime(:,izprime,candidatelocation(Khigh));
    mpolbprimeprime = interp1(vgridomega',mpolbprimetemp,...
        mpolomegaprime,"linear","extrap");
    mpolbprimeprime(mpolbprimeprime<pblimit) = pblimit;

    mprime = (mpolomegaprime - mpolkprimeprime - qprime*mpolbprimeprime)/(wprime*zprime); 
    nprime = (-peta*mprime + sqrt((peta*mprime).^2+4*peta))/(2*peta);

    cprime = wprime.*zprime*nprime + mpolomegaprime - mpolkprimeprime - qprime*mpolbprimeprime;
    cprime(cprime<=0) = 1e-10;

    muprime = 1./cprime;                
    mexpectation =  mexpectation + (1+rprime).*muprime.*repmat(mtransz(:,izprime)',pnumgridomega,1)*mtransA(iA,iAprime);
    mexpectation2 =  mexpectation2 + muprime.*repmat(mtransz(:,izprime)',pnumgridomega,1)*mtransA(iA,iAprime);
    
    else

    K2Lprime = Kprime/tsupplyL(ifuture);
    rprime   = palpha*Aprime*(K2Lprime)^(palpha-1)-pdelta;
    wprime   = (1-palpha)*Aprime*(K2Lprime)^(palpha);        
    qprime   = tq(ifuture);
    
    mpolomegaprime = mpolkprime(:,:,itrans)*(1+rprime) + mpolbprime(:,:,itrans);

    mpolkprimeprime = interp1(vgridomega',squeeze(mpolkprime(:,izprime,ifuture)),...
            mpolomegaprime,"linear","extrap");
    mpolkprimeprime(mpolkprimeprime<0) = 0;
    mpolbprimeprime = interp1(vgridomega',squeeze(mpolbprime(:,izprime,ifuture)),...
            mpolomegaprime,"linear","extrap");
    mpolbprimeprime(mpolbprimeprime<pblimit) = pblimit;

    mprime = (mpolomegaprime - mpolkprimeprime - qprime*mpolbprimeprime)/(wprime*zprime); 
    nprime = (-peta*mprime + sqrt((peta*mprime).^2+4*peta))/(2*peta);

    cprime = wprime.*zprime*nprime + mpolomegaprime - mpolkprimeprime - qprime*mpolbprimeprime;
    cprime(cprime<=0) = 1e-10;

    muprime = 1./cprime;                
    mexpectation =  mexpectation + (1+rprime).*muprime.*repmat(mtransz(:,izprime)',pnumgridomega,1)*mtransA(iA,futureShock);
    mexpectation2 =  mexpectation2 + muprime.*repmat(mtransz(:,izprime)',pnumgridomega,1)*mtransA(iA,futureShock);
        
    end
    end
end

mexpectation = pbeta*mexpectation;
mexpectation2 = pbeta*mexpectation2;

c = 1./(mexpectation+ mlambda(:,:,itrans));
n = (w*mgridz./(peta*c)).^pfrisch;

mlambda_newtemp = 1./mpolc(:,:,itrans) - mexpectation;
mlambda_newtemp(mlambda_newtemp<0) = 0;

mphi_newtemp = q./mpolc(:,:,itrans) - mexpectation2;
mphi_newtemp(mphi_newtemp<0)=0;

mpolkprime_newtemp = w.*mgridz.*n + vgridomega' - mpolc(:,:,itrans) - q*mpolbprime(:,:,itrans);
mpolbprime_newtemp = (w.*mgridz.*n + vgridomega' - c - mpolkprime_newtemp)/q;

mpolkprime_newtemp(mpolkprime_newtemp<=0) = 0;
mpolkprime_newtemp(mlambda(:,:,itrans)>0) = 0;
mlambda_newtemp(mpolkprime_newtemp>0) = 0;

mpolbprime_newtemp(mpolbprime_newtemp<=pblimit) = pblimit;
mpolbprime_newtemp(mphi(:,:,itrans)>0) = pblimit;
mphi_newtemp(mpolbprime_newtemp>pblimit) = 0;


% update
mpolkprime_new(:,:,itrans) = mpolkprime_newtemp;
mpolbprime_new(:,:,itrans) = mpolbprime_newtemp;
mlambda_new(:,:,itrans) = mlambda_newtemp;
mphi_new(:,:,itrans) = mphi_newtemp;
mpolc_new(:,:,itrans) = c;
mpoln(:,:,itrans) = n;

end

% update consumption policy
mpolc = mpolc_new;

% for conservative update
% mpolc = mpolc * weightold + mpolc_new * (1-weightold);

%=========================    
% 3. non-stochastic simulation
%=========================       
% initial distribution
currentdist = startingdist;
tK_new = zeros(size(tsimpath));
tq_new = zeros(size(tsimpath));
tlambda_new = zeros(size(tsimpath));
tphi_new = zeros(size(tsimpath));
tsupplyL_new = zeros(size(tsimpath));
trisky = zeros(size(tsimpath));
triskylow = zeros(size(tsimpath));
triskyhigh = zeros(size(tsimpath));
tklow = zeros(size(tsimpath));
tkhigh = zeros(size(tsimpath));
tBlow = zeros(size(tsimpath));
tBhigh = zeros(size(tsimpath));
riskyratio = zeros(size(mpolkprime));
tK_new(1) = ss.endoK;
tempSupplyL = zeros(size(currentdist));
tr(pathlength+1) = ss.r;

for itrans = 1:pathlength
    
iA = tsimpath(itrans);
nextdist = zeros(size(currentdist));
rprime = tr(itrans+1);
mpolomegaprime = mpolkprime_new(:,:,itrans)*(1+rprime) + mpolbprime_new(:,:,itrans);

for iz = 1:pnumgridz
    
    tempSupplyL(:,iz) = squeeze(mgridz(:,iz).*mpoln(:,iz,itrans));
    
for iomega = 1:pnumgridomega

    nextomega = mpolomegaprime(iomega,iz);
    
    lb = sum(vgridomega<nextomega);
    lb(lb<=0) = 1;
    lb(lb>=pnumgridomega) = pnumgridomega-1;
    ub = lb+1;
    weightlb = (vgridomega(ub) - nextomega)/(vgridomega(ub)-vgridomega(lb));
    weightlb(weightlb<0) = 0;
    weightlb(weightlb>1) = 1;                    
    weightub = 1-weightlb;

    mass = currentdist(iomega,iz);

    for izprime = 1:pnumgridz
        
        nextdist(lb,izprime) = ...
            nextdist(lb,izprime)...
            +mass*mtransz(iz,izprime)*weightlb;

        nextdist(ub,izprime) = ...
            nextdist(ub,izprime)...
            +mass*mtransz(iz,izprime)*weightub;
    end 
end
end

%update

tOmega(itrans) = sum(vgridomega'.*sum(currentdist,2));       
tB(itrans+1) = sum(mpolbprime(:,:,itrans).*currentdist,'all');
tK_new(itrans+1) = sum(mpolkprime(:,:,itrans).*currentdist,'all');
tsupplyL_new(itrans)= sum(tempSupplyL.*currentdist,'all');        
tY(itrans) = vgridA(tsimpath(itrans))*tK_new(itrans)^palpha*tsupplyL(itrans)^(1-palpha);
tC(itrans) = sum(squeeze(mpolc(:,:,itrans)).*currentdist,'all');   
tlambda(itrans) = sum(mlambda(:,:,itrans).*currentdist,'all');
tlambda_new(itrans) = sum(mlambda_new(:,:,itrans).*currentdist,'all');
tphi(itrans) = sum(mphi(:,:,itrans).*currentdist,'all');
tphi_new(itrans) = sum(mphi_new(:,:,itrans).*currentdist,'all');

A = vgridA(iA);
w = (1-palpha)*A*(tK_new(itrans)/tsupplyL_new(itrans))^(palpha);

% riskyratio(:,:,itrans+1) = (mpolkprime(:,:,itrans).*(1+tr(itrans+1))./(mpolkprime(:,:,itrans).*(1+tr(itrans+1))+mpolbprime(:,:,itrans)));
% trisky(itrans+1) = sum(riskyratio(:,:,itrans+1).*currentdist,'all');
riskyratio(:,:,itrans+1) = (mpolkprime(:,:,itrans)./(mpolkprime(:,:,itrans)+mpolbprime(:,:,itrans)));
trisky(itrans+1) = sum(riskyratio(:,:,itrans+1).*currentdist,'all');
currentdistlow = currentdist;
currentdistlow(:,3:5) = 0;
currentdisthigh = currentdist;
currentdisthigh(:,1:3) = 0;
triskylow(itrans+1) = sum(riskyratio(:,:,itrans+1).*currentdistlow,'all')./sum(currentdistlow,'all');
triskyhigh(itrans+1) = sum(riskyratio(:,:,itrans+1).*currentdisthigh,'all')./sum(currentdisthigh,'all');
tklow(itrans+1) = sum(mpolkprime(:,:,itrans).*currentdistlow,'all')./sum(currentdistlow,'all');
tBlow(itrans+1) = sum(mpolbprime(:,:,itrans).*currentdistlow,'all')./sum(currentdistlow,'all');
tkhigh(itrans+1) = sum(mpolkprime(:,:,itrans).*currentdisthigh,'all')./sum(currentdisthigh,'all');
tBhigh(itrans+1) = sum(mpolbprime(:,:,itrans).*currentdisthigh,'all')./sum(currentdisthigh,'all');
tq_new(itrans) = (sum( currentdist.*( mgridomega + mpoln(:,:,itrans) .* mgridz * w ...
               - mpolc_new(:,:,itrans) - mpolkprime_new(:,:,itrans) ...
               + pBbar*tq(itrans)),'all'))/pBbar;

% alternatively,
% tq_new(itrans) = (tB(itrans+1)+pBbar*tq(itrans))/pBbar;

currentdist = nextdist;

end
    
%=========================  
% check the convergence and update the price
%=========================  

% clearing (convergence) conditions

error2  = mean(abs(...
        [tK(burnin+1:pathlength-burnin)-tK_new(burnin+1:pathlength-burnin) ;...
        tsupplyL(burnin+1:pathlength-burnin)-tsupplyL_new(burnin+1:pathlength-burnin);...
        tq(burnin+1:pathlength-burnin)-tq_new(burnin+1:pathlength-burnin) ;...
        tlambda(burnin+1:pathlength-burnin)-tlambda_new(burnin+1:pathlength-burnin); ...
        tphi(burnin+1:pathlength-burnin)-tphi_new(burnin+1:pathlength-burnin)...
        ].^2));

errorK      = tK - tK_new;

tK          = weightold1.*tK            +(1-weightold1).*tK_new;
mlambda     = weightold2.*mlambda       +(1-weightold2).*mlambda_new;
mphi        = weightold2.*mphi          +(1-weightold2).*mphi_new;
mpolkprime  = weightold3.*mpolkprime    +(1-weightold3).*mpolkprime_new;
mpolbprime  = weightold3.*mpolbprime    +(1-weightold3).*mpolbprime_new;
tsupplyL    = weightold4.*tsupplyL      +(1-weightold4).*tsupplyL_new;
tq          = weightold5.*tq            +(1-weightold5).*tq_new;

% if pnumiter_ge == 1
% error2 = 1;
% end

if verbose2 == true && (floor((pnumiter_ge-1)/20) == (pnumiter_ge-1)/20) || error2<= tol_ge
%=========================  
% interim report
%=========================  
fprintf(' \n');
fprintf('Market Clearing Results \n');
fprintf('Error: %.10f \n', error2);

% plot
subplot(2,2,1);
plot(1:pathlength,tK(1:pathlength));hold on;
plot(1:pathlength,tK_new(1:pathlength),'-.');
xlim([1,pathlength]);
hold off;
legend("Predicted","Realized","location","northeast");

subplot(2,2,2);
plot(1:pathlength,tsupplyL(1:pathlength));hold on;
plot(1:pathlength,tsupplyL_new(1:pathlength),'-.');
xlim([1,pathlength]);
hold off;
legend("Predicted","Realized","location","northeast");

subplot(2,2,3);
plot(1:pathlength,tq(1:pathlength));hold on;
plot(1:pathlength,tq_new(1:pathlength),'-.');
xlim([1,pathlength]);
hold off;
legend("Predicted","Realized","location","northeast");

subplot(2,2,4);
plot(1:pathlength,tB(1:pathlength));hold on;
xlim([1,pathlength]);
ylim([-1,1]);
hold off;
legend("Excess bond demand","location","northeast");

pause(0.1)         
fprintf('\n');
toc;

%=========================  
% save (mid)
%=========================  
save '../solutions/WIP_ks1997endolabor_bc.mat';
end

pnumiter_ge = pnumiter_ge+1;
% error2(pnumiter_ge==2) = 1e-30;
end

%=========================  
% save (final)
%=========================  
save '../solutions/ks1997endolabor_bc.mat';


%%
%=========================  
% final report
%========================= 
% define investment
tI = tK([2:pathlength,pathlength]) - tK([1,1:(pathlength-1)])*(1-pdelta);

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
fprintf('mean log(output): %.4f \n', mean(log(tY)));
fprintf('st. dev. log(output): %.4f \n', std(log(tY)));
fprintf('skewness log(output): %.4f \n', skewness(log(tY)));
fprintf('------------------------ \n');
fprintf('mean log(investment): %.4f \n', mean(log(tI)));
fprintf('st. dev. log(investment): %.4f \n', std(log(tI)));
fprintf('skewness log(investment): %.4f \n', skewness(log(tI)));
fprintf('------------------------ \n');
fprintf('mean log(consumption): %.4f \n', mean(log(tC)));
fprintf('st. dev. log(consumption): %.4f \n', std(log(tC)));
fprintf('skewness log(consumption): %.4f \n', skewness(log(tC)));

fprintf('\n');
fprintf('======================== \n');
fprintf('Business cycle statistics for the HP-filtered time series\n');
fprintf('======================== \n');
[~,vYhpfilter] = hpfilter(log(tY),1600);
[~,vIhpfilter] = hpfilter(log(tI),1600);
[~,vChpfilter] = hpfilter(log(tC),1600);
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
endostate = tK(burnin+1:end-2);
exostate = vgridA(tsimpath(burnin+1:end-1))';
endostateprime = tK(burnin+2:end-1);
intercept = ones(size(endostate));

% independent variable
x = [intercept ...
    ,log(endostate) ...
    ,log(exostate) ...
    ,log(endostate).*log(exostate)...
    ];

% dependent variable
y = log(endostateprime);

[coeff,bint1,r1,rint1,R1] = regress(y,x);
fprintf('======================== \n');
fprintf('Fitting the true LoM into the log-linear specification\n');
fprintf('======================== \n');
disp(['R-squared: ',num2str(R1(1))]);
fprintf(' \n');
fitlm(x,y,'Intercept',false)

% recover the implied dynamics
startingpoint = burnin+1;
startingendo = endostate(1);
recovered = ones(1,pathlength - burnin)*startingendo;
for itrans = 1:(pathlength - burnin-1)
% endoStateTemp = endoState(iTrans);
endostatetemp = recovered(itrans);
exostatetemp = exostate(itrans);
tempX = [1 ...
    ,log(endostatetemp)...
    ,log(exostatetemp)...
    ,log(endostatetemp)*log(exostatetemp)...
    ];
tempVal = coeff'*tempX';
recovered(itrans+1) = exp(tempVal);
end

samplePeriod = 500:1000;
figure;
plot(samplePeriod,endostate(samplePeriod),'Color','red','LineWidth',1.5);hold on;
plot(samplePeriod,recovered(samplePeriod),'Color','blue','LineStyle','--','LineWidth',1.5);
legend("True LoM","Linear LoM","location","best","FontSize",15)
location = ['../figures/lom.pdf'];
saveas(gcf, location);


%%
%=========================  
% fitting wage dynamics into the linear specification 
%========================= 
tw   = (1-palpha)*vgridA(tsimpath)'.*(tK(1:end-1)./tsupplyL).^(palpha);

endostate = tK(burnin+1:end-2);
exostate = vgridA(tsimpath(burnin+1:end-1))';
endoprice = tw(burnin+1:end-1);
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
recovered = zeros(1,pathlength - burnin);
for itrans = 1:(pathlength - burnin-1)
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



%%
%=========================  
% fitting bond price dynamics into the linear specification 
%========================= 
endostate = tK(burnin+1:end-2);
exostate = vgridA(tsimpath(burnin+1:end-1))';
endoprice = tq(burnin+1:end-1);
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
recovered = zeros(1,pathlength - burnin-1);
for itrans = 1:(pathlength - burnin-1)
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
location = ['../figures/lom_q.pdf'];
saveas(gcf, location);

% save for q figures
save '../solutions/ks1997endolabor_bc.mat';
