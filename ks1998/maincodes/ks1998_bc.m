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
load '../solutions/ks1998_ss.mat';
ss = load('../solutions/ks1998_ss.mat');

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
vgridz      = [0.25,1.00];
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
T = 1000;
burnin = 100;
% T = 2001;% T = 1001; % T = 3001;% T = 5001;% T = 10001;
% burnin = 500;
pathlength = T+burnin;
pinitialpoint = 1;
tsimpath = fnSimulator(pinitialpoint,mtransA,burnin+T);    

%=========================
%labor supply (exogenous)
%=========================
[vL,~] = eigs(mtransz',1);
vL = vL/sum(vL);
disp(vL);
supplyL = zeros(1,2);
supplyL(1) = 2*(vgridz*vL(1:2));
supplyL(2) = 2*(vgridz*vL(3:4));

%=========================        
% declare equilibrium objects
%=========================    
mpolc = zeros(pnumgrida,pnumgridz,pathlength);
mpolc_new = zeros(pnumgrida,pnumgridz,pathlength);
mlambda   = zeros(pnumgrida,pnumgridz,pathlength);
mlambda_new  = zeros(pnumgrida,pnumgridz,pathlength);
mpolaprime = zeros(pnumgrida,pnumgridz,pathlength);
mPolaprime_new = zeros(pnumgrida,pnumgridz,pathlength);

%=========================            
% start and end points
%=========================    
startingdist = currentdist;
for itrans = 1:pathlength
mpolaprime(:,:,itrans) = ss.mpolaprime;
mlambda(:,:,itrans) = ss.mlambda;
mpolc(:,:,itrans) = ss.mpolc;
end
    
%=========================     
% initial guess
%========================= 
tK = ss.K*ones(pathlength+1,1)+ normrnd(0,0.00001,pathlength+1,1);
% for the initial iteration, slightly perturbed capital path is necessary.
tY = zeros(size(pathlength,1));
tC = zeros(size(pathlength,1));
tlambda =  sum(ss.mlambda.*ss.currentdist,"all")*ones(pathlength,1);

%%
%=========================    
% resume from the last one
%=========================    
% use the following line if you wish to start from where you stopped
% before.
% load '../solutions/WIP_ks1998_bc.mat';

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
K2L = K/supplyL(iA);                
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

    K2Lprime = Kprime/supplyL(iAprime);
    rprime   = palpha*Aprime*(K2Lprime)^(palpha-1)-pdelta;
    wprime   = (1-palpha)*Aprime*(K2Lprime)^(palpha);    

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
    
    mpolaprime_temp = weightlow*mpolaprime(:,izprime,candidatelocation(Klow)) ...
                + (1-weightlow)*mpolaprime(:,izprime,candidatelocation(Khigh));
    mpolaprimeprime = interp1(vgrida',mpolaprime_temp,...
        squeeze(mpolaprime(:,:,itrans)),"linear","extrap");
    cprime = wprime.*zprime + (1+rprime).*squeeze(mpolaprime(:,:,itrans)) - mpolaprimeprime;
    cprime(cprime<=1e-10) = 1e-10;
    muprime = 1./cprime;                
    mexpectation =  mexpectation + (1+rprime).*muprime.*repmat(temptrans(:,izprime)',pnumgrida,1)*mtransA(iA,iAprime);
    
    else

    % for the realized future shock level on the simulated path
    mpolaprimeprime = interp1(vgrida',squeeze(mpolaprime(:,izprime,ifuture)),...
            squeeze(mpolaprime(:,:,itrans)),"linear","extrap");
    cprime = wprime.*zprime + (1+rprime).*squeeze(mpolaprime(:,:,itrans)) - mpolaprimeprime;
    cprime(cprime<=1e-10) = 1e-10;
    muprime = 1./cprime;                
    mexpectation =  mexpectation + (1+rprime).*muprime.*repmat(temptrans(:,izprime)',pnumgrida,1)*mtransA(iA,futureShock);
        
    end
    
end

mlambda_newtemp = 1./mpolc(:,:,itrans) - pbeta*mexpectation;
mexpectation = pbeta*mexpectation + mlambda(:,:,itrans);
c = 1./mexpectation;
mpolaprime_newtemp = w.*mgridz + (1+r).*mgrida - c;
mlambda_newtemp(mpolaprime_newtemp>vgridamin) = 0;
c = (c - (vgridamin-mpolaprime_newtemp)).*(mpolaprime_newtemp<=vgridamin) + c.*(mpolaprime_newtemp>vgridamin);
mpolaprime_newtemp(mpolaprime_newtemp<=vgridamin) = vgridamin;

% update
mpolaprime_new(:,:,itrans) = mpolaprime_newtemp;
mlambda_new(:,:,itrans) = mlambda_newtemp;
mpolc_new(:,:,itrans) = c;

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
tK_new(1) = sum(vgrida*sum(currentdist,2));

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
tY(itrans) = vgridA(tsimpath(itrans))*tK_new(itrans)^palpha*supplyL(tsimpath(itrans))^(1-palpha);
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
        tlambda(burnin+1:pathlength-burnin)-tlambda_new(burnin+1:pathlength-burnin)].^2));

errorK      = tK - tK_new;
tK          = weightold1.*tK            +(1-weightold1).*tK_new;
mlambda     = weightold2.*mlambda       +(1-weightold2  ).*mlambda_new;
mpolaprime  = weightold3.*mpolaprime    +(1-weightold3).*mpolaprime_new;
 
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
plot(1:pathlength,tK(1:pathlength));hold on;
plot(1:pathlength,tK_new(1:pathlength),'-.');
xlim([1,pathlength]);
hold off;
legend("Predicted","Realized","location","northeast");

pause(0.1)         
fprintf('\n');
toc;

%=========================  
% save (mid)
%=========================  
save '../solutions/WIP_ks1998_bc.mat';
end

pnumiter_ge = pnumiter_ge+1;
end

%=========================  
% save (final)
%=========================  
save '../solutions/ks1998_bc.mat';


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
recovered = ones(1,pathlength- burnin)*startingendo;
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

samplePeriod = 500:800;
figure;
plot(samplePeriod,endostate(samplePeriod),'Color','red','LineWidth',1.5);hold on;
plot(samplePeriod,recovered(samplePeriod),'Color','blue','LineStyle','--','LineWidth',1.5);
legend("True LoM","Linear LoM","location","best","FontSize",15);
ylabel("Aggregate capital stock","FontSize",15);
xlabel("Time (quarter)","FontSize",15);
ax = gca;
box off;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.Grid off;
location = ['../figures/lom.pdf'];
saveas(gcf, location);

%%
% plot - for report
predictedpath = tK(burnin+1:end-2);
realizedpath = tK_new(burnin+1:end-2);

samplePeriod = 500:800;
figure;
plot(samplePeriod,predictedpath(samplePeriod),'LineWidth',2.0);hold on;
plot(samplePeriod,realizedpath(samplePeriod),'-.','LineWidth',2.0);
plot(samplePeriod,recovered(samplePeriod),':','Color',"black",'LineWidth',2.0);
xlim([min(samplePeriod),max(samplePeriod)]);
hold off;
box off;
legend("Predicted","Realized","Krusell and Smith (1998)","location","northeast");
ylabel("Aggregate capital stock","FontSize",15);
xlabel("Time (quarter)","FontSize",15);
ax = gca;
ax.FontSize = 15; 
set(gcf, 'PaperPosition', [0 0 6 5]); %Position plot at left hand corner with width a and height b.
set(gcf, 'PaperSize', [6 5]); %Set the paper to have width a and height b.Grid off;
location = ['../figures/rtm_kscomp.pdf'];
saveas(gcf, location);
