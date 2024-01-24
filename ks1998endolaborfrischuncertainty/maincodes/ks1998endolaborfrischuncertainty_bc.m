%% Krusell and Smith (1998) with endogenous labor supply and aggregate uncertainty shock
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
load '../solutions/ks1998endolaboruncertainty_ss.mat';
ss = load('../solutions/ks1998endolaboruncertainty_ss.mat');

%=========================
% numerical parameters
%=========================
tol_ge          = 1e-6;
weightold1      = 0.9500;
weightold2      = 0.9500;
weightold3      = 0.9500;
weightold4      = 0.9500;

%=========================
% aggregate shock
%=========================
% aggregate productivity shock
pnumgridA   = 2;
mtransA     = [0.875,0.125;0.125,0.875];
vgridAL     = [0.9916,1.0084];
vgridAH     = -(1-vgridAL)*1.6+[1.00,1.00]; % Bloom et al. (2018)

% transition matrix from Bloom et al. (2018)
pnumgridS = 2;
mtransS = [1-0.026,0.026;1-0.943,0.943];
vgridS  = [0,1];

% stationary distribution of the transition shock
[vS,~] = eigs(mtransS',1);
vS = vS/sum(vS);
disp(vS);
disp(vS(1)*vgridAL+vS(2)*vgridAH); % average grid

% combine the aggregate shocks
vgridS          = kron(vgridS,ones(pnumgridA,1));
vgridSlocation  = kron((1:1:pnumgridS)',ones(pnumgridA,1));
vgridAlocation  = kron(ones(pnumgridS,1),(1:1:pnumgridA)');
mtransA_orig    = mtransA;
mtransA         = kron(mtransS,mtransA_orig);
pnumgridA       = pnumgridA*pnumgridS;
vgridAextended  = [vgridAL,vgridAH];

%=========================
% simulation path
%=========================
T = 2001;% T = 1001;% T = 2001;% T = 5001;% T = 10001;
burnin = 500;
pathlength = T+burnin;
pinitialpoint = 1;

% simulation
seed = 100;
rng(seed);
tsimpathA = fnSimulator(pinitialpoint,mtransA_orig,burnin+T);

seed = 1000;
rng(seed);
tsimpathS = fnSimulator(pinitialpoint,mtransS,burnin+T);

tsimpath = (tsimpathA==1).*(tsimpathS==1)*1 ...
         + (tsimpathA==2).*(tsimpathS==1)*2 ...
         + (tsimpathA==1).*(tsimpathS==2)*3 ...
         + (tsimpathA==2).*(tsimpathS==2)*4;

% alternative
% seed = 100;
% rng(seed);
% tsimpath    = fnSimulator(pinitialpoint,mtransA,burnin+T);

%=========================        
% declare equilibrium objects
%=========================    
mpolc       = zeros(pnumgrida,pnumgridz,pathlength);
mpolaprime  = zeros(pnumgrida,pnumgridz,pathlength);
mlambda     = zeros(pnumgrida,pnumgridz,pathlength);
mlambda_new = zeros(pnumgrida,pnumgridz,pathlength);
mpolc_new   = zeros(pnumgrida,pnumgridz,pathlength);
mpoln       = zeros(size(mpolc));
mpolaprime_new = zeros(pnumgrida,pnumgridz,pathlength);

%=========================            
% start and end points
%=========================    
startingdist = currentdist;
for itrans = 1:pathlength
mpolaprime(:,:,itrans) = ss.mpolaprime;
mpoln(:,:,itrans) = ss.mpoln;
mlambda(:,:,itrans) = ss.mlambda;
mpolc(:,:,itrans) = ss.mpolc;
end

%=========================     
% initial guess
%========================= 
tK = ss.K*ones(pathlength+1,1)+ normrnd(0,0.0001,pathlength+1,1);
% for the initial iteration, slightly perturbed capital path is necessary.
tsupplyL = ss.supplyL*ones(pathlength,1);
tY = zeros(pathlength,1);
tC = zeros(pathlength,1);
tlambda =  sum(ss.mlambda.*ss.currentdist,"all")*ones(pathlength,1);

%%
%=========================    
% resume from the last one
%=========================    
% use the following line if you wish to start from where you stopped
% before.
% load '../solutions/WIP_ks1998endolaboruncertainty_bc.mat';

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

iAll = tsimpath(itrans);
ipast = itrans-1;
ipast(ipast<1) = 1;
iAllpast = tsimpath(ipast);
iA = vgridAlocation(iAll);
iS = vgridSlocation(iAll);
iSpast = vgridSlocation(iAllpast);
vgridA = vgridAL*(iSpast==1) + vgridAH*(iSpast==2);
vgridAfuture = vgridAL*(iS==1) + vgridAH*(iS==2);

A = vgridA(iA);
K = tK(itrans);
supplyL = tsupplyL(itrans);

%given K, all the prices are known.
r   = palpha*A*(K/supplyL)^(palpha-1)-pdelta;
mmu = r+pdelta;
w   = (1-palpha)*A*(K/supplyL)^(palpha);

if itrans == pathlength
iAllfuture = tsimpath(itrans);
ifuture = itrans;
Kprime = tK(itrans+1);
else
iAllfuture = tsimpath(itrans+1);
ifuture = itrans+1;
Kprime = tK(itrans+1);
end

% expected future value (rational expectation)
mpolaprimecurrent = squeeze(mpolaprime(:,:,itrans));
mexpectation = 0;
for iAllprime = 1:pnumgridA
    
    iAprime = vgridAlocation(iAllprime);
    Aprime = vgridAfuture(iAprime);

    for izprime = 1:pnumgridz
    zprime = vgridz(izprime);
    
    if iAllfuture ~= iAllprime || itrans == pathlength
    % find a period where the future shock realization is the same as
    % iAprime and the capital stock is closest to Kprime from below and above.
    candidate = tK(find(tsimpath==iAllprime)); % the candidates with the specific exogenous state
    candidatelocation = find(tsimpath==iAllprime); % the candidates' location
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
    
    mpolaprimetemp = weightlow*mpolaprime(:,izprime,candidatelocation(Klow)) ...
                + (1-weightlow)*mpolaprime(:,izprime,candidatelocation(Khigh));
    mpolaprimeprime = interp1(vgrida',mpolaprimetemp,...
                mpolaprimecurrent,"linear","extrap");

    K2Lprimelow = Kprime/tsupplyL(candidatelocation(Klow));
    rprimelow   = palpha*Aprime*(K2Lprimelow)^(palpha-1)-pdelta;
    wprimelow   = (1-palpha)*Aprime*(K2Lprimelow)^(palpha);        

    K2Lprimehigh = Kprime/tsupplyL(candidatelocation(Khigh));
    rprimehigh   = palpha*Aprime*(K2Lprimehigh)^(palpha-1)-pdelta;
    wprimehigh   = (1-palpha)*Aprime*(K2Lprimehigh)^(palpha);        

    rprime      = weightlow*rprimelow + (1-weightlow)*rprimehigh;
    wprime      = weightlow*wprimelow + (1-weightlow)*wprimehigh;

    mprime = ((1+rprime)*mpolaprime(:,:,itrans) - mpolaprimeprime)/(wprime*zprime); 
    nprime = (-peta*mprime + sqrt((peta*mprime).^2+4*peta))/(2*peta);

    cprime = wprime.*zprime*nprime + (1+rprime).*mpolaprimecurrent - mpolaprimeprime;
    cprime(cprime<=0) = 1e-10;

    muprime = 1./cprime;                
    mexpectation =  mexpectation + (1+rprime).*muprime.*repmat(mtransz(:,izprime)',pnumgrida,1)*mtransA(iAll,iAllprime);
    
    else

    mpolaprimeprime = interp1(vgrida',squeeze(mpolaprime(:,izprime,ifuture)),...
                mpolaprimecurrent,"linear","extrap");
    K2Lprime = Kprime/tsupplyL(ifuture);
    rprime   = palpha*Aprime*(K2Lprime)^(palpha-1)-pdelta;
    wprime   = (1-palpha)*Aprime*(K2Lprime)^(palpha);        

    mprime = ((1+rprime)*mpolaprime(:,:,itrans) - mpolaprimeprime)/(wprime*zprime); 
    nprime = (-peta*mprime + sqrt((peta*mprime).^2+4*peta))/(2*peta);

    cprime = wprime.*zprime*nprime + (1+rprime).*mpolaprimecurrent - mpolaprimeprime;
    cprime(cprime<=0) = 1e-10;

    muprime = 1./cprime;                
    mexpectation =  mexpectation + (1+rprime).*muprime.*repmat(mtransz(:,izprime)',pnumgrida,1)*mtransA(iAll,iAllfuture);
        
    end
    end
end

mexpectation = pbeta*mexpectation;
c = 1./(mexpectation+ mlambda(:,:,itrans));
n = (w*mgridz./(peta*c)).^pfrisch;
mlambda_newtemp = 1./mpolc(:,:,itrans) - mexpectation;
mpolaprime_newtemp = w.*mgridz.*n + (1+r).*mgrida - c;
mlambda_newtemp(mpolaprime_newtemp>vgridamin) = 0;
c = (c - (vgridamin-mpolaprime_newtemp)).*(mpolaprime_newtemp<=vgridamin) + c.*(mpolaprime_newtemp>vgridamin);
mpolaprime_newtemp(mpolaprime_newtemp<=vgridamin) = vgridamin;

% update
mpolaprime_new(:,:,itrans) = mpolaprime_newtemp;
mlambda_new(:,:,itrans) = mlambda_newtemp;
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
tlambda_new = zeros(size(tsimpath));
tsupplyL_new = zeros(size(tsimpath));
tK_new(1) = sum(vgrida*sum(currentdist,2));
tempSupplyL = zeros(size(currentdist));

for itrans = 1:pathlength
    
iAll = tsimpath(itrans);
iA = vgridAlocation(iAll);
iS = vgridSlocation(iAll);
vgridA = vgridAL*(iS==1) + vgridAH*(iS==2);

nextdist = zeros(size(currentdist));

for iz = 1:pnumgridz
    
    tempSupplyL(:,iz) = squeeze(mgridz(:,iz).*mpoln(:,iz,itrans));
    
for ia = 1:pnumgrida

    a = vgrida(ia);
    nexta = mpolaprime_new(ia,iz,itrans);
    
    lb = sum(vgrida<nexta);
    lb(lb<=0) = 1;
    lb(lb>=pnumgrida) = pnumgrida-1;
    ub = lb+1;
    weightlb = (vgrida(ub) - nexta)/(vgrida(ub)-vgrida(lb));
    weightlb(weightlb<0) = 0;
    weightlb(weightlb>1) = 1;                    
    weightub = 1-weightlb;

    mass = currentdist(ia,iz);

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
tK_new(itrans+1) = sum(vgrida'.*sum(nextdist,2));        
tsupplyL_new(itrans)= sum(tempSupplyL.*currentdist,'all');        
tY(itrans) = vgridA(iA)*tK_new(itrans)^palpha*tsupplyL(itrans)^(1-palpha);
tC(itrans) = sum(squeeze(mpolc(:,:,itrans)).*currentdist,'all');   
tlambda(itrans) = sum(mlambda(:,:,itrans).*currentdist,'all');
tlambda_new(itrans) = sum(mlambda_new(:,:,itrans).*currentdist,'all');
currentdist = nextdist;

end
    
%=========================  
% check the convergence and update the price
%=========================  
% market clearing
error2  = mean(abs(...
        [tK(burnin+1:pathlength-burnin)-tK_new(burnin+1:pathlength-burnin) ;...
        tsupplyL(burnin+1:pathlength-burnin)-tsupplyL_new(burnin+1:pathlength-burnin);...
        tlambda(burnin+1:pathlength-burnin)-tlambda_new(burnin+1:pathlength-burnin)].^2));

errorK      = tK - tK_new;
tK          = weightold1.*tK            +(1-weightold1).*tK_new;
mlambda     = weightold2.*mlambda       +(1-weightold2).*mlambda_new;
mpolaprime  = weightold3.*mpolaprime    +(1-weightold3).*mpolaprime_new;
tsupplyL    = weightold4.*tsupplyL      +(1-weightold4).*tsupplyL_new;
 
% if pnumiter_ge == 1
% error2 = 1;
% end

if verbose2 == true && (floor((pnumiter_ge-1)/10) == (pnumiter_ge-1)/10) || error2<= tol_ge
%=========================  
% interim report
%=========================  
fprintf(' \n');
fprintf('Market Clearing Results \n');
fprintf('Error: %.10f \n', error2);

% plot
subplot(1,2,1);
plot(1:pathlength,tK(1:pathlength));hold on;
plot(1:pathlength,tK_new(1:pathlength),'-.');
xlim([1,pathlength]);
hold off;
legend("Predicted","Realized","location","northeast");

subplot(1,2,2);
plot(1:pathlength,tsupplyL(1:pathlength));hold on;
plot(1:pathlength,tsupplyL_new(1:pathlength),'-.');
xlim([1,pathlength]);
hold off;
legend("Predicted","Realized","location","northeast");


pause(0.1)         
fprintf('\n');
toc;

%=========================  
% save (mid)
%=========================  
save '../solutions/WIP_ks1998endolaboruncertainty_bc.mat';
end

pnumiter_ge = pnumiter_ge+1;
end

%=========================  
% save (final)
%=========================  
save '../solutions/ks1998endolaboruncertainty_bc.mat';


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
endostate = tK(burnin+1:end-1);
exostate = vgridAextended(tsimpath(burnin+1:end-1))';
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
tw   = (1-palpha)*vgridAextended(tsimpath)'.*(tK./tsupplyL).^(palpha);

endostate = tK(burnin+1:end-2);
exostate = vgridAextended(tsimpath(burnin+1:end-1))';
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
