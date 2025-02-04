%% Krusell and Smith (1997) with endogenous labor supply
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
ss = load('../solutions/ks1997endolabor_ss.mat');
globalSol = load('../solutions/ks1997endolabor_bc.mat');
load('../solutions/ks1997endolabor_bc.mat');

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
mpolc_new(:,:,itrans) = c;

end

EE = (mpolc_new - mpolc)./mpolc;

%%
% the individual household state needs to be specified.
% the monotonicity result is robust over the choice of different individual states.

iz = floor(pnumgridz/2);
ik = floor(pnumgridomega/2);

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