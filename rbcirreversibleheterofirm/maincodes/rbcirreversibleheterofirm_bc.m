%% An RBC model with heterogeneous firms and irreversible investments
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
partialeq       = false;        % turning on/off ge effect

%=========================
% numerical parameters
%=========================
tol_ge          = 1e-7;
weightold1      = 0.9500;
weightold2      = 0.9980;

%=========================
% load the stationary equilibrium allocations
%=========================
load '../solutions/rbcirreversibleheterofirm_ss.mat';
ss = load('../solutions/rbcirreversibleheterofirm_ss.mat');

% exogenous boundaries just in case of initial explosive path
pl              = ss.eq.p*0.1;%1.5;
pu              = ss.eq.p*4.0;%3.5;
Kl              = ss.eq.k*0.1;%0.5;
Ku              = ss.eq.k*7.0;%5.0;

%=========================
% aggregate shock
%=========================
% parameter setup
prho_A = 0.859;
psigma_A = 0.014;

% Tauchen method
pnumgridA = 5;
[mtransA, vgridA] = ...
fnTauchen(prho_A, 0, psigma_A^2, pnumgridA, 2);
vgridA = exp(vgridA);

%=========================
% simulation path
%=========================
seed = 100;
rng(seed);
T = 1001;
% T = 2001;% T = 3001;% T = 5001;% T = 10001;
burnin = 500;
pathlength = T+burnin;
pinitialpoint = 1;
tsimpath = fnSimulator(pinitialpoint,mtransA,burnin+T);    

%=========================        
% declare equilibrium objects
%=========================    
tjderiv     = zeros(pnumgridk,pnumgridz,pathlength);
tj1deriv    = zeros(pnumgridk,pnumgridz,pathlength);
tpoleffdvd  = zeros(pnumgridk,pnumgridz,pathlength);
tpoli       = zeros(pnumgridk,pnumgridz,pathlength);
tpoll       = zeros(pnumgridk,pnumgridz,pathlength);
tpoly       = zeros(pnumgridk,pnumgridz,pathlength);
tpolk       = zeros(pnumgridk,pnumgridz,pathlength);
tlambda     = zeros(pnumgridk,pnumgridz,pathlength);

%=========================  
% compute the equilibrium allocations
%=========================  
eq.p = zeros(pathlength,1);
eq.w = zeros(pathlength,1);
eq.k = zeros(pathlength,1);
eq.c = zeros(pathlength,1);

%=========================            
% start and end points
%=========================    
startingdist = currentdist;
endvalue = ss.mjderiv;
for it = 1:pathlength
    tjderiv(:,:,it) = endvalue;
    tlambda(:,:,it) = ss.mlambda;
end
    
%=========================     
% initial guess
%========================= 
tp = ss.eq.p*ones(pathlength,1);
tK = ss.eq.k*ones(pathlength+1,1)+ normrnd(0,0.000001,pathlength+1,1);
% for the initial iteration, slightly perturbed capital path is necessary.
      
%%
%=========================    
% resume from the last one
%=========================    
% use the following line if you wish to start from where you stopped
% before.
% load '../solutions/rbcirreversibleheterofirm_bc.mat';

%%        
%=========================    
% outer loop
%=========================    
tic;
error2 = 10;
pnumiter_ge = 1;

while error2>tol_ge
   
% macro setup: wage
tw = peta./tp;
  
%=========================    
% 2. backward solution
%=========================    
for itrans = pathlength:-1:1
        
    iA = tsimpath(itrans);
    A = vgridA(iA);            
    w = tw(itrans);
    p = tp(itrans);
    K = tK(itrans);
    
    if itrans == pathlength
    Kprime = tK(itrans+1);
    futureShock = tsimpath(itrans);
    ifuture = itrans;
    else
    Kprime = tK(itrans+1);
    futureShock = tsimpath(itrans+1);
    ifuture = itrans+1;
    end
    
    % expected future value (rational expectation)
    tempjderiv = 0;
    % risk-free interest rate
    eq.r(itrans) = 0;
    for iAprime = 1:pnumgridA
    if futureShock ~= iAprime || itrans == pathlength
    % find a period where the future shock realization is the same as
    % iAprime and the capital stock is closest to Kprime from below and above.
    candidate = tK(find(tsimpath==iAprime)); % the candidates with the specific exogenous state
    candidatelocation = find(tsimpath==iAprime); % the candidates' location
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
    
    tempjderiv = tempjderiv + mtransA(iA,iAprime)*...
                    (weightlow*tjderiv(:,:,candidatelocation(Klow)) ...
                  + (1-weightlow)*tjderiv(:,:,candidatelocation(Khigh)));
    
    eq.r(itrans) = eq.r(itrans)+mtransA(iA,iAprime)*...
                   (weightlow*tp(candidatelocation(Klow)) ...
                 + (1-weightlow)*tp(candidatelocation(Khigh)));
    
    else 

    % for the realized future shock level on the simulated path
    tempjderiv = tempjderiv + mtransA(iA,futureShock)*tjderiv(:,:,ifuture);
    eq.r(itrans) = eq.r(itrans) + mtransA(iA,futureShock)*tp(ifuture);        
    
    end
    end
    
    tj0deriv = tempjderiv*mtransz';
    eq.r(itrans) = (1/pbeta)*(p)/eq.r(itrans) - 1;
    z = mgridz.*A;
    k = mgridk;                 
    l = ((pgamma)*z.*k.^(palpha)/w).^(1/(1-pgamma));
    
    for iz = 1:pnumgridz
        
        tempmat = tj0deriv(:,iz)';
        kprimeopt = fnoptinvirreversible(pbeta,vgridk,p,tempmat,tlambda(:,iz,itrans));
        kprimeopt = (threshold + (1-pdelta)*vgridk ).*(kprimeopt<threshold+(1-pdelta)*vgridk) ...
                  + (kprimeopt                     ).*(kprimeopt>=threshold+(1-pdelta)*vgridk);

        low = sum(mgridkk <= kprimeopt',2);
        low(low<=0) = 1;
        low(low>=pnumgridk) = pnumgridk-1;
        high = low + 1;

        indxlow     = sub2ind(size(mgridkk),(1:pnumgridk)',low);
        indxhigh    = sub2ind(size(mgridkk),(1:pnumgridk)',high);
        
        weightdown  = (mgridkk(indxhigh) - kprimeopt')...
                    ./(mgridkk(indxhigh)-mgridkk(indxlow));
        weightup    = 1-weightdown;

        tempfuturevalderiv = tj0deriv(low ,iz).*weightdown + ...
                             tj0deriv(high,iz).*weightup;
        lambda =  p - pbeta*tempfuturevalderiv;
        lambda(lambda<0) = 0;
        
        mkprimeopt(:,iz) = kprimeopt;
        tlambda(:,iz,itrans) = lambda;
        %mlambda(:,iz) = 0;

    end 

    % updated value
    optvalderiv     = p*(1-pgamma)*(pgamma/w).^(pgamma/(1-pgamma)).*z.^(1/(1-pgamma))...
                    .* (palpha/(1-pgamma)).*k.^(palpha/(1-pgamma)-1) ...
                    + (1-pdelta).*(p-tlambda(:,:,itrans));
    
    % update
    tj1deriv(:,:,itrans)  = optvalderiv;
    tpoleffdvd(:,:,itrans) = z.*k.^(palpha).*l.^(pgamma) - (mkprimeopt-(1-pdelta)*k);
    tpoli(:,:,itrans)     = mkprimeopt - (1-pdelta)*k;
    tpoll(:,:,itrans)     = l;
    tpoly(:,:,itrans)     = z.*k.^(palpha).*l.^(pgamma);
    tpolk(:,:,itrans)     = mkprimeopt;                
    
    if verbose1 == true && (floor((itrans-1)/30) == (itrans-1)/30)
    dispthis=['Iteration is in process : ', num2str(itrans),' from ',num2str(pathlength)];
    disp(dispthis)
    toc;   
    end
    
end

tjderiv = tj1deriv; 

% for conservative update
% tj = tj1 * weightold + tj0 * (1-weightold);

%=========================    
% 3. non-stochastic simulation
%=========================       
% initial distribution
currentdist = startingdist;

for itrans = 1:pathlength
% empty distribution to save the next distribution
nextdist = zeros(size(currentdist)); 

for iz = 1:pnumgridz
for ik = 1:pnumgridk

nextk = tpolk(ik,iz,itrans);
lb = sum(vgridk<nextk);
lb(lb<=0) = 1;
lb(lb>=pnumgridk) = pnumgridk-1;
ub = lb+1;
weightlb = (vgridk(ub) - nextk)./(vgridk(ub)-vgridk(lb));
weightlb(weightlb<0) = 0;
weightlb(weightlb>1) = 1;                    
weightub = 1-weightlb;

mass = currentdist(ik,iz);

for futureiz = 1:pnumgridz
    
    nextdist(lb,futureiz) = ...
        nextdist(lb,futureiz)...
        +mass.*mtransz(iz,futureiz).*weightlb;

    nextdist(ub,futureiz) = ...
        nextdist(ub,futureiz)...
        +mass.*mtransz(iz,futureiz).*weightub;

end                 
end
end

%=========================  
% 4. compute the equilibrium allocations
%=========================  
eq.p(itrans) = tp(itrans);
eq.w(itrans) = tw(itrans); 

eq.k(itrans) = sum(repmat(vgridk',1,pnumgridz).*currentdist,'all');
eq.i(itrans) = sum(tpoli(:,:,itrans).*currentdist,'all');
eq.ik(itrans) = sum(tpoli(:,:,itrans)./repmat(vgridk',1,pnumgridz).*currentdist,'all');
eq.l(itrans) = sum(tpoll(:,:,itrans).*currentdist,'all');
eq.y(itrans) = sum(tpoly(:,:,itrans).*currentdist,'all');
tempc        = sum(tpoleffdvd(:,:,itrans).*currentdist,'all');
eq.c(itrans) = tempc*(tempc>1e-8) + 1e-8*(tempc<1e-8);

currentdist = nextdist;
    
end

% for the last period's future capital
eq.k(itrans+1) = sum(repmat(vgridk',1,pnumgridz).*currentdist,'all');

%=========================  
% market clearing 
%=========================  
eq.p_implied = 1./eq.c;
if partialeq==true
    eq.p_implied = ones(size(eq.p_implied))./ss.eq.c;
end

% exogenous boundary
eq.p_implied = (eq.p_implied<pl)*pl+(eq.p_implied>pu)*pu+(eq.p_implied>=pl).*(eq.p_implied<=pu).*eq.p_implied;
eq.k = (eq.k<Kl)*Kl+(eq.k>Ku)*Ku+(eq.k>=Kl).*(eq.k<=Ku).*eq.k;

% market clearing
errorp = abs(eq.p_implied - eq.p);        
errorK = abs(eq.k - tK);

%=========================  
% 5. update
%=========================  
oldGridp = tp;
tK = tK*weightold1+eq.k*(1-weightold1);
tp = tp*weightold2+eq.p_implied*(1-weightold2);
 
if pnumiter_ge == 1
error2 = 1;
else
error2 = mean([errorp;errorK].^2);
end

if verbose2 == true && (floor((pnumiter_ge-1)/10) == (pnumiter_ge-1)/10) || error2<= tol_ge
%=========================  
% interim report
%=========================  
fprintf(' \n');
fprintf('Market Clearing Results \n');
fprintf('Error: %.10f \n', error2);
fprintf('Max error: %.10f \n', max([errorp;errorK]));
fprintf('Mean squared error: %.10f \n', mean([errorp;errorK].^2)); 

% plot
subplot(1,2,1);
plot(1:pathlength,tp(1:pathlength));hold on;
plot(1:pathlength,eq.p_implied(1:pathlength),'-.');
xlim([1,pathlength]);
hold off;
legend("Predicted p","Realized p","location","northeast");

subplot(1,2,2);
plot(1:pathlength,tK(1:pathlength));hold on;
plot(1:pathlength,eq.k(1:pathlength),'-.');
xlim([1,pathlength]);
hold off;
legend("Predicted K","Realized K","location","northeast");

pause(0.1)         
fprintf('\n');
toc;

%=========================  
% save (mid)
%=========================  
save '../solutions/WIP_rbcirreversibleheterofirm_bc.mat';

end

pnumiter_ge = pnumiter_ge+1;

end

%=========================  
% save (final)
%=========================  
save '../solutions/rbcirreversibleheterofirm_bc.mat';


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
fprintf('mean log(output): %.4f \n', mean(log(eq.y)));
fprintf('st. dev. log(output): %.4f \n', std(log(eq.y)));
fprintf('skewness log(output): %.4f \n', skewness(log(eq.y)));
fprintf('------------------------ \n');
fprintf('mean log(investment): %.4f \n', mean(log(eq.i)));
fprintf('st. dev. log(investment): %.4f \n', std(log(eq.i)));
fprintf('skewness log(investment): %.4f \n', skewness(log(eq.i)));
fprintf('------------------------ \n');
fprintf('mean log(consumption): %.4f \n', mean(log(eq.c)));
fprintf('st. dev. log(consumption): %.4f \n', std(log(eq.c)));
fprintf('skewness log(consumption): %.4f \n', skewness(log(eq.c)));

fprintf('\n');
fprintf('======================== \n');
fprintf('Business cycle statistics for the HP-filtered time series\n');
fprintf('======================== \n');
[~,vYhpfilter] = hpfilter(log(eq.y),1600);
[~,vIhpfilter] = hpfilter(log(eq.i),1600);
[~,vChpfilter] = hpfilter(log(eq.c),1600);
fprintf('mean log(output): %.4f \n', mean((vYhpfilter)));
fprintf('st. dev. log(output): %.4f \n', std((vYhpfilter)));
fprintf('skewness log(output): %.4f \n', skewness((vYhpfilter)));
fprintf('------------------------ \n');
fprintf('mean log(investment): %.4f \n', mean((vIhpfilter)));
fprintf('st. dev. log(investment): %.4f \n', std((vIhpfilter)));
fprintf('skewness log(investment): %.4f \n', skewness((vIhpfilter)));
fprintf('------------------------ \n');
fprintf('mean log(consumption): %.4f \n', mean((vChpfilter)));
fprintf('st. dev. log(consumption): %.4f \n', std((vChpfilter)));
fprintf('skewness log(consumption): %.4f \n', skewness((vChpfilter)));

%%
%=========================  
% dynamic consistency report
%========================= 
fprintf('\n');
fprintf('======================== \n');
fprintf('dynamic consistency report for rtm \n');
fprintf('======================== \n');
disp(['max absolute error (in pct. of steady-state K): ', num2str(100*max(abs(errorK))/ss.eq.k),'%']);
disp(['root mean sqaured error (in pct. of steady-state K): ', num2str(100*sqrt(mean(errorK.^2))/ss.eq.k),'%']);
fprintf('\n');
figure;
hist(100*errorK/ss.eq.k,100);
xlim([-1,1]);
xlabel("Dynamic consistency error (in pct. of steady-state K)")
location = ['../figures/err_hist.pdf'];
saveas(gcf, location);

%%
%=========================  
% fitting LoM into the linear specification 
%========================= 
endostate = tK(burnin+1:end-2);
exostate = vgridA(tsimpath(burnin+1:end-1));
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
recovered = ones(1,pathlength- burnin)*startingendo;
for iTrans = 1:(pathlength - burnin-1)
% endoStateTemp = endoState(iTrans);
endostatetemp = recovered(iTrans);
exostatetemp = exostate(iTrans);
tempX = [1 ...
    ,log(endostatetemp)...
    ,log(exostatetemp)...
    ,log(endostatetemp)*log(exostatetemp)...
    ];
tempVal = coeff'*tempX';
recovered(iTrans+1) = exp(tempVal);
end

samplePeriod = 500:1000;
figure;
plot(samplePeriod,endostate(samplePeriod),'Color','red','LineWidth',1.5);hold on;
plot(samplePeriod,recovered(samplePeriod),'Color','blue','LineStyle','--','LineWidth',1.5);
legend("True LoM","Linear LoM","location","best","FontSize",15)
location = ['../figures/lom.pdf'];
saveas(gcf, location);


%%
% plot - for report
predictedpath = log(tK(burnin+1:end-2));
realizedpath = log(eq.k(burnin+1:end-2));

samplePeriod = 500:1000;
figure;
plot(samplePeriod,predictedpath(samplePeriod),'LineWidth',2.0);hold on;
plot(samplePeriod,realizedpath(samplePeriod),'-.','LineWidth',2.0);
plot(samplePeriod,log(recovered(samplePeriod)),':','Color',"black",'LineWidth',2.0);
xlim([min(samplePeriod),max(samplePeriod)]);
hold off;
box off;
legend("Predicted","Realized","Linear LoM","location","northeast");
ylabel("log(K_{t})","FontSize",15);
xlabel("Time (quarter)","FontSize",15);
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.
location = ['../figures/rtm_comp.pdf'];
saveas(gcf, location);

%%
%=========================  
% fitting price dynamics into the linear specification 
%========================= 
endostate = tK(burnin+1:end-2);
exostate = vgridA(tsimpath(burnin+1:end-1));
endoprice = tp(burnin+1:end-1);
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
xlim([min(samplePeriod),max(samplePeriod)]);
legend("True LoM","Linear LoM","location","best","FontSize",15);
hold off;
box off;
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.
location = ['../figures/lom_p.pdf'];
saveas(gcf, location);
