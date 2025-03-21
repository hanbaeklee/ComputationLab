%% Solving the model in Krusell and Smith (1998)
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
load '../solutions/ks1998endolabor_ss.mat';
ss = load('../solutions/ks1998endolabor_ss.mat');

%=========================
% numerical parameters
%=========================
tol_ge          = 1e-6;
weightold1      = 0.9000;
weightold2      = 0.9000;
weightold3      = 0.9000;

%=========================
% aggregate shock
%=========================
% idiosyncratic income shock
pnumgridz   = 2;
vgridz      = [0.00,1.00];
mtransz     = [0.525,0.350,0.03125,0.09375;...
               0.035,0.84,0.0025,0.1225;...
               0.09375,0.03125,0.292,0.583;...
               0.0099,0.1151,0.0245,0.8505];      

% aggregate productivity shock
pnumgridA   = 2;
mtransA = [sum(sum(mtransz(1:2,1:2))),sum(sum(mtransz(1:2,3:4))); ...
           sum(sum(mtransz(3:4,1:2))),sum(sum(mtransz(3:4,3:4)))];
mtransA = mtransA/pnumgridA;
vgridA = [0.99,1.01];

%=========================
% simulation path
%=========================
seed = 100;
rng(seed);
T = 2001;% T = 1001; % T = 3001;% T = 5001;% T = 10001;
burnin = 500;
pathlength = T+burnin;
pinitialpoint = 1;
tsimpath = fnSimulator(pinitialpoint,mtransA,burnin+T);    

%=========================        
% declare equilibrium objects
%=========================    
mpolc = zeros(pnumgrida,pnumgridz,pathlength);
mpolc_new = zeros(pnumgrida,pnumgridz,pathlength);
mpoln = zeros(pnumgrida,pnumgridz,pathlength);
mpoln_new = zeros(pnumgrida,pnumgridz,pathlength);
mlambda   = zeros(pnumgrida,pnumgridz,pathlength);
mlambda_new  = zeros(pnumgrida,pnumgridz,pathlength);
mpolaprime = zeros(pnumgrida,pnumgridz,pathlength);
mpolaprime_new = zeros(pnumgrida,pnumgridz,pathlength);

%=========================            
% start and end points
%=========================    
startingdist = currentdist;
for itrans = 1:pathlength
mpolaprime(:,:,itrans) = ss.mpolaprime;
mlambda(:,:,itrans) = ss.mlambda;
mpolc(:,:,itrans) = ss.mpolc;
mpoln(:,:,itrans) = ss.mpoln;
end
    
%=========================     
% initial guess
%========================= 
tK = ss.K*ones(pathlength+1,1)+ normrnd(0,0.000001,pathlength+1,1);
% for the initial iteration, slightly perturbed capital path is necessary.
tsupplyL = ss.supplyL*ones(pathlength,1);
tY = zeros(size(pathlength,1));
tC = zeros(size(pathlength,1));
tlambda = sum(ss.mlambda.*ss.currentdist,"all")*ones(pathlength,1);

%%
%=========================    
% resume from the last one
%=========================    
% use the following line if you wish to start from where you stopped
% before.
load '../solutions/WIP_ks1998endolabor_bc.mat';

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
K2L = K/supplyL;                
range = 1:2;
if iA == 2
    range = 3:4;
end

r   = palpha*A*(K2L)^(palpha-1)-pdelta;
w   = (1-palpha)*A*(K2L)^(palpha);

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
for ishockprime = 1:pnumgridz*pnumgridA
    
    iAprime = (ishockprime<=2)*1 + (ishockprime>=3)*2;
    Aprime = vgridA(iAprime);
    izprime = ishockprime;
    rangeprime = 1:2;
    if iAprime == 2
    izprime = izprime - pnumgridz;
    rangeprime = 3:4;
    end
    temptrans = (mtransz(range,rangeprime)./sum(mtransz(range,rangeprime),2));
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

    mpolaprime_temp = weightlow*mpolaprime(:,izprime,candidatelocation(Klow)) ...
                + (1-weightlow)*mpolaprime(:,izprime,candidatelocation(Khigh));
    mpolaprimeprime = interp1(vgrida',mpolaprime_temp,...
            squeeze(mpolaprime(:,:,itrans)),"linear","extrap");
    
    nprime_temp = weightlow*mpoln(:,izprime,candidatelocation(Klow)) ...
           + (1-weightlow)*mpoln(:,izprime,candidatelocation(Khigh));
    nprime = interp1(vgrida',nprime_temp,...
            squeeze(mpolaprime(:,:,itrans)),"linear","extrap");
    
    cprime = wprime.*zprime.*nprime + (1+rprime).*squeeze(mpolaprime(:,:,itrans)) - mpolaprimeprime;
    cprime(cprime<=1e-10) = 1e-10;
    muprime = 1./cprime;                
    mexpectation =  mexpectation + (1+rprime).*muprime.*repmat(temptrans(:,izprime)',pnumgrida,1)*mtransA(iA,iAprime);
    
    else

    % for the realized future shock level on the simulated path
    mpolaprimeprime = interp1(vgrida',squeeze(mpolaprime(:,izprime,ifuture)),...
             squeeze(mpolaprime(:,:,itrans)),"linear","extrap");
    nprime   = interp1(vgrida',squeeze(mpoln(:,izprime,ifuture)),...
             squeeze(mpolaprime(:,:,itrans)),"linear","extrap");
    K2Lprime = Kprime/tsupplyL(ifuture);
    rprime   = palpha*Aprime*(K2Lprime)^(palpha-1)-pdelta;
    wprime   = (1-palpha)*Aprime*(K2Lprime)^(palpha);        
    cprime   = wprime.*zprime.*nprime + (1+rprime).*squeeze(mpolaprime(:,:,itrans)) - mpolaprimeprime;
    cprime(cprime<=1e-10) = 1e-10;
    muprime = 1./cprime;                
    mexpectation =  mexpectation + (1+rprime).*muprime.*repmat(temptrans(:,izprime)',pnumgrida,1)*mtransA(iA,futureShock);
        
    end
    
end

mexpectation = pbeta*mexpectation;
c = 1./(mexpectation+ mlambda(:,:,itrans));
n = 1-(1-ptheta).*mpolc(:,:,itrans)./(w.*mgridz.*ptheta);
% mlambda_newtemp = 1./mpolc(:,:,itrans) - pbeta*mexpectation;
mlambda_newtemp = 1./(w.*mgridz.*mpoln(:,:,itrans) + (1+r).*mgrida - mpolaprime(:,:,itrans)) - pbeta*mexpectation;
mlambda_newtemp(mlambda_newtemp<0) = 0;
mpolaprime_newtemp = w.*mgridz.*mpoln(:,:,itrans) + (1+r).*mgrida - c;
mlambda_newtemp(mpolaprime_newtemp>vgridamin) = 0;
c = (c - (vgridamin-mpolaprime_newtemp)).*(mpolaprime_newtemp<=vgridamin) + c.*(mpolaprime_newtemp>vgridamin);
mpolaprime_newtemp(mpolaprime_newtemp<=vgridamin) = vgridamin;

% update
mpolaprime_new(:,:,itrans) = mpolaprime_newtemp;
mlambda_new(:,:,itrans) = mlambda_newtemp;
mpolc_new(:,:,itrans) = c;
mpoln_new(:,:,itrans) = n;
mpoln_new(mpoln_new<0) = 0;
mpoln_new(mpoln_new>1) = 1;

end

% update consumption policy
mpolc = mpolc_new;
% for conservative update
% mpolc = mpolc * weightold + mpolc_new * (1-weightold);
% mpoln = mpoln * weightold + mpoln_new * (1-weightold);

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
    
iA = tsimpath(itrans);
range = 1:2;
if iA == 2
    range = 3:4;
end
ifuture = itrans+1;
ifuture(ifuture>pathlength)=pathlength;
iAprime = tsimpath(ifuture);
rangeprime = 1:2;
if iAprime == 2
rangeprime = 3:4;
end
temptrans = (mtransz(range,rangeprime)./sum(mtransz(range,rangeprime),2));

%basic setup
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
    weightUB = 1-weightlb;

    mass = currentdist(ia,iz);

    for izprime = 1:pnumgridz
        
        nextdist(lb,izprime) = ...
            nextdist(lb,izprime)...
            +mass*temptrans(iz,izprime)*weightlb;

        nextdist(ub,izprime) = ...
            nextdist(ub,izprime)...
            +mass*temptrans(iz,izprime)*weightUB;
    end 
end
end

%update
tK_new(itrans+1) = sum(vgrida'.*sum(nextdist,2));        
tsupplyL_new(itrans)= sum(tempSupplyL.*currentdist,'all');        
tY(itrans) = vgridA(tsimpath(itrans))*tK_new(itrans)^palpha*tsupplyL(itrans)^(1-palpha);
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
        tsupplyL(burnin+1:pathlength-burnin)-tsupplyL_new(burnin+1:pathlength-burnin)].^2));

errorK      = tK - tK_new;
tK          = weightold1.*tK            +(1-weightold1).*tK_new;
tsupplyL    = weightold2.*tsupplyL      +(1-weightold2).*tsupplyL_new;
mlambda     = weightold3.*mlambda       +(1-weightold3).*mlambda_new;
mpolaprime  = weightold4.*mpolaprime    +(1-weightold4).*mpolaprime_new;
mpoln       = weightold5.*mpoln         +(1-weightold5).*mpoln_new;

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
save '../solutions/WIP_ks1998endolabor_bc.mat';
end

pnumiter_ge = pnumiter_ge+1;
end

%=========================  
% save (final)
%=========================  
save '../solutions/ks1998endolabor_bc.mat';


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
hold off;
box off;
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.Grid off;
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
hold off;
box off;
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.Grid off;
location = ['../figures/lom_w.pdf'];
saveas(gcf, location);

