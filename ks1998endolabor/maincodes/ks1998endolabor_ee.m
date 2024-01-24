%% Solving the model in Krusell and Smith (1998)
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
ss = load('../solutions/ks1998_ss.mat');
globalSol = load('../solutions/ks1998_bc.mat');
load('../solutions/ks1998_bc.mat');

%%
%=========================
%backward solution
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
    Kprime = K;

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

    if futureShock ~= iAprime        
        
    candidate = tK(find(tsimpath==iAprime));
    candidateLocation = find(tsimpath==iAprime);
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
    
    mpolaprime_temp = weightlow*mpolaprime(:,izprime,candidateLocation(Klow)) ...
                + (1-weightlow)*mpolaprime(:,izprime,candidateLocation(Khigh));
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