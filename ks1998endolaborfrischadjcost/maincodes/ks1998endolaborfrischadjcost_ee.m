%% Krusell and Smith (1998) with endogenous labor supply and convex adjustment cost
% 2023.09.25
% Hanbaek Lee (hanbaeklee1@gmail.com)
% When you use the code, please cite the paper 
% "A Dynamically Consistent Global Nonlinear Solution 
% Method in the Sequence Space and Applications."
%=========================    
% this file is to compute the euler equation errors.
%=========================    
%=========================    
% housekeeping
%=========================
% clc;
% clear variables;
% close all; 
% fnPath = './functions';
% addpath(fnPath);

%=========================
%load ss
%=========================
ss = load('../solutions/ks1998endolaborfrischadjcost_ss.mat');
globalSol = load('../solutions/ks1998endolaborfrischadjcost_bc.mat');
load('../solutions/ks1998endolaborfrischadjcost_bc.mat');

%%
%=========================
%backward solution
%=========================
for itrans = pathlength:-1:1

iA = tsimpath(itrans);
A = vgridA(iA);
K = tK(itrans);
supplyL = tsupplyL(itrans);

%given K, all the prices are known.
r   = palpha*A*(K/supplyL)^(palpha-1)-pdelta;
mmu = r+pdelta;
w   = (1-palpha)*A*(K/supplyL)^(palpha);

if itrans == pathlength

    futureShock = tsimpath(itrans);
    iFuture = itrans;
    Kprime = K;

else

    futureShock = tsimpath(itrans+1);
    iFuture = itrans+1;
    Kprime = tK(itrans+1);

end


% expected future value (rational expectation)
mexpectation = 0;
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
    
    mpolaprimetemp = weightlow*mpolaprime(:,izprime,candidatelocation(Klow)) ...
                + (1-weightlow)*mpolaprime(:,izprime,candidatelocation(Khigh));
    mpolaprimeprime = interp1(vgrida',mpolaprimetemp,...
        squeeze(mpolaprime(:,:,itrans)),"linear","extrap");

    K2Lprimelow = Kprime/tsupplyL(candidatelocation(Klow));
    rprimelow   = palpha*Aprime*(K2Lprimelow)^(palpha-1)-pdelta;
    wprimelow   = (1-palpha)*Aprime*(K2Lprimelow)^(palpha);        

    K2Lprimehigh = Kprime/tsupplyL(candidatelocation(Khigh));
    rprimehigh   = palpha*Aprime*(K2Lprimehigh)^(palpha-1)-pdelta;
    wprimehigh   = (1-palpha)*Aprime*(K2Lprimehigh)^(palpha);        
    
    psi2        = -(pmu/2)*((mpolaprimeprime./mpolaprime(:,:,itrans)).^2-1);  
    rprime      = weightlow*rprimelow + (1-weightlow)*rprimehigh;
    wprime      = weightlow*wprimelow + (1-weightlow)*wprimehigh;

    mprime = ((1+rprime)*mpolaprime(:,:,itrans) - (pmu/2).*((mpolaprimeprime-mpolaprime(:,:,itrans))./mpolaprime(:,:,itrans)).^2.*mpolaprime(:,:,itrans)...
           - mpolaprimeprime)/(wprime*zprime); 
    nprime = (-peta*mprime + sqrt((peta*mprime).^2+4*peta))/(2*peta);

    cprime = wprime.*zprime*nprime + (1+rprime).*squeeze(mpolaprime(:,:,itrans)) - mpolaprimeprime ...
           - (pmu/2).*((mpolaprimeprime-mpolaprime(:,:,itrans))./mpolaprime(:,:,itrans)).^2.*mpolaprime(:,:,itrans);
    cprime(cprime<=0) = 1e-10;

    muprime = 1./cprime;                
    mexpectation =  mexpectation + (1+rprime-psi2).*muprime.*repmat(mtransz(:,izprime)',pnumgrida,1)*mtransA(iA,iAprime);
    
    else

    mpolaprimeprime = interp1(vgrida',squeeze(mpolaprime(:,izprime,ifuture)),...
            squeeze(mpolaprime(:,:,itrans)),"linear","extrap");
    K2Lprime = Kprime/tsupplyL(ifuture);
    rprime   = palpha*Aprime*(K2Lprime)^(palpha-1)-pdelta;
    wprime   = (1-palpha)*Aprime*(K2Lprime)^(palpha);     
    psi2     = -(pmu/2)*((mpolaprimeprime./mpolaprime(:,:,itrans)).^2-1);  

    mprime = ((1+rprime)*mpolaprime(:,:,itrans) - (pmu/2).*((mpolaprimeprime-mpolaprime(:,:,itrans))./mpolaprime(:,:,itrans)).^2.*mpolaprime(:,:,itrans)...
           - mpolaprimeprime)/(wprime*zprime); 
    nprime = (-peta*mprime + sqrt((peta*mprime).^2+4*peta))/(2*peta);

    cprime = wprime.*zprime*nprime + (1+rprime).*squeeze(mpolaprime(:,:,itrans)) - mpolaprimeprime ...
           - (pmu/2).*((mpolaprimeprime-mpolaprime(:,:,itrans))./mpolaprime(:,:,itrans)).^2.*mpolaprime(:,:,itrans);
    cprime(cprime<=0) = 1e-10;

    muprime = 1./cprime;                
    mexpectation =  mexpectation + (1+rprime - psi2).*muprime.*repmat(mtransz(:,izprime)',pnumgrida,1)*mtransA(iA,futureShock);
        
    end
    end
end

mexpectation = pbeta*mexpectation;
mexpectation = (mexpectation + mlambda(:,:,itrans))./(1+pmu*(mpolaprime(:,:,itrans)-mgrida)./mgrida); 
c = 1./(mexpectation);
n = (w*mgridz./(peta*c)).^pfrisch;
% mlambdaNewTemp = 1./(w.*mGridz.*n + (1+r).*mgrida - mpolaprime(:,:,iTrans)) - mExpectation;
mlambda_newtemp = 1./mpolc(:,:,itrans) - mexpectation;
mpolaprime_newtemp = w.*mgridz.*n + (1+r).*mgrida - c - (pmu/2).*((mpolaprime(:,:,itrans)-mgrida)./mgrida).^2.*mgrida;
mlambda_newtemp(mpolaprime_newtemp>vgridamin) = 0;
c = (c - (vgridamin-mpolaprime_newtemp)).*(mpolaprime_newtemp<=vgridamin) + c.*(mpolaprime_newtemp>vgridamin);
mpolaprime_newtemp(mpolaprime_newtemp<=vgridamin) = vgridamin;

% update
mpolaprime_new(:,:,itrans) = mpolaprime_newtemp;
mlambda_new(:,:,itrans) = mlambda_newtemp;
mpolc_new(:,:,itrans) = c;
mpoln(:,:,itrans) = n;

end

EE = (mpolc_new - mpolc)./mpolc;

%%
% the individual household state needs to be specified.
% the monotonicity result is robust over the choice of different individual states.

iz = floor(pnumgridz/2);
ik = floor(pnumgrida/2);

ts1 = tK(burnin+1:end-1);%(vSimPath==3);
ts2 = squeeze(EE(ik,iz,burnin+1:end)); %(vSimPath==3);

figure;
scatter(1:length(ts2),log(abs(ts2)));
xlabel("Time (year)","FontSize",15);
ylabel("Euler equation error (in Log)","FontSize",15);
hold off;
location = ['../figures/ee_timeseries.pdf'];
saveas(gcf, location);

figure;
hist(log(abs(ts2)),100);
xlabel("Euler equation error (in Log)","FontSize",15);
ylabel("Distribution","FontSize",15);
location = ['../figures/ee_hist.pdf'];
saveas(gcf, location);

figure;
scatter(ts1,log(abs(ts2)));
xlabel("Aggregate capital stock","FontSize",15);
ylabel("Euler equation error (in Log)","FontSize",15);
location = ['../figures/ee.pdf'];
saveas(gcf, location);