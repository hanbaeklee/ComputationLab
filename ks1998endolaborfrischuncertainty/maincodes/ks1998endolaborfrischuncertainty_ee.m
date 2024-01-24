%% Krusell and Smith (1998) with endogenous labor supply
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
ss = load('../solutions/ks1998endolaboruncertainty_ss.mat');
globalSol = load('../solutions/ks1998endolaboruncertainty_bc.mat');
load('../solutions/ks1998endolaboruncertainty_bc.mat');

%%
%=========================
%backward solution
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
Kprime = K;
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
    
    if iAllfuture ~= iAllprime        
    candidate = tK(find(tsimpath==iAllprime));
    candidateLocation = find(tsimpath==iAllprime);
    candidate(candidateLocation<burnin) = [];
    candidateLocation(candidateLocation<burnin) = [];
    [candidate,index] = sort(candidate);
    candidateLocation = candidateLocation(index);

    Klow = sum(candidate<Kprime);
    Klow(Klow<=1) = 1;
    Klow(Klow>=length(index)) = length(index)-1;
    Khigh = Klow+1;
    weightlow = (candidate(Khigh) - Kprime)/(candidate(Khigh)-candidate(Klow));
    weightlow(weightlow<0) = 0;
    weightlow(weightlow>1) = 1;
    
    mpolaprimetemp = weightlow*mpolaprime(:,izprime,candidateLocation(Klow)) ...
                + (1-weightlow)*mpolaprime(:,izprime,candidateLocation(Khigh));
    mpolaprimeprime = interp1(vgrida',mpolaprimetemp,...
                mpolaprimecurrent,"linear","extrap");

    K2Lprimelow = Kprime/tsupplyL(candidateLocation(Klow));
    rprimelow   = palpha*Aprime*(K2Lprimelow)^(palpha-1)-pdelta;
    wprimelow   = (1-palpha)*Aprime*(K2Lprimelow)^(palpha);        

    K2Lprimehigh = Kprime/tsupplyL(candidateLocation(Khigh));
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

EE = (mpolc_new - mpolc)./mpolc;

%%
% the individual household state needs to be specified.
% the monotonicity result is robust over the choice of different individual states.

iz = floor(pnumgridz/2);
ik = floor(pnumgrida/2);

ts1 = tK;%(vSimPath==3);
ts2 = squeeze(EE(ik,iz,:)); %(vSimPath==3);

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