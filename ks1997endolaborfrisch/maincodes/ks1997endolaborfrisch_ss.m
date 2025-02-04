%% Krusell and Smith (1997) with endogenous labor supply
% 2023.09.25
% Hanbaek Lee (hanbaeklee1@gmail.com)
% When you use the code, please cite the paper 
% "A Dynamically Consistent Global Nonlinear Solution 
% Method in the Sequence Space and Applications."
%=========================    
% this file is to compute the stationary equilibrium.
%=========================    
%=========================    
% housekeeping
%=========================
clc;
clear variables;
close all; 
fnPath = '../functions';
addpath(fnPath);

%=========================
% macro choices
%=========================
eigvectormethod  = 2;            % eigenvector method on/off for stationary dist computation
verbose          = true;         % ge loop interim reports on/off

%=========================
% parameters
%=========================
palpha      = 0.36; 
pbeta       = 0.99; 
pdelta      = 0.025;
prho        = 0.90;
psigma      = 0.05;
pfrisch     = 1.00;
peta        = 8.00;
pblimit     = -2.40;
pBbar       = 1.00;

%=========================
% numerical parameters - grids
%=========================
%idiosyncratic income shock
pnumgridz   = 5;
[mtransz, vgridz] = ...
fnTauchen(prho, 0, psigma^2, pnumgridz, 2);
vgridz = exp(vgridz');            
if pnumgridz ==1
    mtransz = 1;
end

%wealth grid
pnumgridomega   = 200;

%finer grid near smaller wealth (Maliar, Maliar, and Valli, 2010)
vgridomegamin   = pblimit;
vgridomegamax   = 300;
x               = linspace(0,0.5,pnumgridomega);
y               = x.^5/max(x.^5);
vgridomega      = vgridomegamin+(vgridomegamax-vgridomegamin)*y;

%=========================
% numerical parameters
%=========================
tol_labor       = 1e-10;
tol_ge          = 1e-10;
tol_pfi         = 1e-8;
tol_hist        = 1e-10;
weightold1      = 0.9900;
weightold2      = 0.9900;
weightold3      = 0.9999;
weightold4      = 0.9500;
weightold5      = 0.9500;

%=========================    
% extra setup
%=========================    
% grids for matrix operations
mgridomega      = repmat(vgridomega',1,pnumgridz);
mgridz          = repmat(vgridz,pnumgridomega,1);
mgridomegaa     = repmat(vgridomega,pnumgridomega,1);

%=========================    
% initial guess
%=========================    
K           = 17.5;
supplyL     = 0.33;
q           = 0.989;%pbeta;%*0.99;
currentdist = ones(pnumgridomega,pnumgridz)/(pnumgridomega*pnumgridz);

%=========================    
% declare equilibrium objects
%=========================
mpolc       = repmat(0.01*vgridomega',1,pnumgridz);
mpolomegaprime = mgridomega;
mpolkprime  = zeros(size(mpolc));
mpolkprime_new = zeros(size(mpolc));
mpolbprime  = zeros(size(mpolc));
mpolbprime_new = zeros(size(mpolc));
mlambda     = zeros(size(mpolc));
mlambda_new = zeros(size(mpolc));
mphi        = zeros(size(mpolc));
mphi_new    = zeros(size(mpolc));


%%
%=========================    
% resume from the last one
%=========================    
% use the following line if you wish to start from where you stopped
% before.
load '../solutions/WIP_ks1997endolabor_ss.mat';

%%
%=========================    
% outer loop
%=========================    
tic;
error2 = 10;
pnumiter_ge = 1;

while error2>tol_ge

% macro setup: given K, all the prices are known.
r   = palpha*(K/supplyL)^(palpha-1)-pdelta;
mmu = r+pdelta;
w   = (1-palpha)*(K/supplyL)^(palpha);

%=========================
% policy function iteration
%=========================
% note that for the stationary equilibrium, I do not run extra loop for the
% policy function iteration. the policy function is simultaneously updated
% with the price (aggregate capital) to boost the speed.

error = 10;
pnumiter_pfi = 1;
while error>tol_pfi
% while pnumiter_pfi<500

mexpectation = 0;
mexpectation2 = 0;
for izprime = 1:pnumgridz
    
    zprime = vgridz(izprime);
    
    rprime = r;
    wprime = w;
    qprime = q;

    mpolkprimeprime = interp1(vgridomega',squeeze(mpolkprime(:,izprime)),...
                    squeeze(mpolomegaprime),"linear","extrap"); 
    mpolbprimeprime = interp1(vgridomega',squeeze(mpolbprime(:,izprime)),...
                    squeeze(mpolomegaprime),"linear","extrap"); 
    mpolkprimeprime(mpolkprimeprime<0) = 0;
    mpolbprimeprime(mpolbprimeprime<pblimit) = pblimit;
    
    mprime = (mpolomegaprime - mpolkprimeprime - qprime*mpolbprimeprime)/(wprime*zprime); 
    nprime = (-peta*mprime + sqrt((peta*mprime).^2+4*peta))/(2*peta);
    
    cprime = wprime.*zprime*nprime + (  mpolomegaprime - mpolkprimeprime - qprime*mpolbprimeprime);
    cprime(cprime<=0) = 1e-10;
    
    muprime = 1./cprime;                
    mexpectation =  mexpectation + repmat(mtransz(:,izprime)',pnumgridomega,1).*(1+rprime).*muprime;
    mexpectation2 =  mexpectation2 + repmat(mtransz(:,izprime)',pnumgridomega,1).*muprime;
    
end

mexpectation = pbeta*mexpectation;
mexpectation2 = pbeta*mexpectation2;

c = 1./(mexpectation+mlambda);
mpoln = (w.*mgridz./(peta*c)).^pfrisch;

mlambda_new = 1./(w.*mgridz.*mpoln + mgridomega - mpolkprime - q*mpolbprime) - mexpectation;
mlambda_new(mlambda_new<0)=0;

mphi_new = q./(w.*mgridz.*mpoln + mgridomega - mpolkprime - q*mpolbprime) - mexpectation2;
mphi_new(mphi_new<0)=0;

mpolkprime_new = w.*mgridz.*mpoln + mgridomega - mpolc - q*mpolbprime;
mpolbprime_new = (w.*mgridz.*mpoln + mgridomega - c - mpolkprime_new)/q;

mpolkprime_new(mpolkprime_new<=0) = 0;
mpolkprime_new(mlambda>0) = 0;
mlambda_new(mpolkprime_new>0) = 0;

mpolbprime_new(mpolbprime_new<=pblimit) = pblimit;
mpolbprime_new(mphi>0) = pblimit;
mphi_new(mpolbprime_new>pblimit) = 0;

mpolc_new = c;
mpolomegaprime_new = (1+r)*mpolkprime_new + mpolbprime_new;

% iteration routine
% error = max([abs(mpolkprime-mpolkprime_new),abs(mpolbprime-mpolbprime_new)].^2,[],"all");
error = max(abs(mpolomegaprime-mpolomegaprime_new).^2,[],"all");

pnumiter_pfi = pnumiter_pfi+1;

mlambda     = weightold4.*mlambda       + (1-weightold4).*mlambda_new;
mphi        = weightold4.*mphi          + (1-weightold4).*mphi_new;
mpolkprime  = weightold5.*mpolkprime    + (1-weightold5).*mpolkprime_new;
mpolbprime  = weightold5.*mpolbprime    + (1-weightold5).*mpolbprime_new;
mpolc       = weightold5.*mpolc         + (1-weightold5).*mpolc_new;
mpolomegaprime = weightold5.*mpolomegaprime    + (1-weightold5).*mpolomegaprime_new;


end

%%
%=========================
if eigvectormethod == 1
%=========================
%option1: histogram-evolving method (non-stochastic)
%=========================
error = 10;
while error>tol_hist
% empty distribution to save the next distribution  
nextdist = zeros(size(currentdist));

for iz = 1:pnumgridz
for ia = 1:pnumgridomega

    a = vgridomega(ia);
    nexta = mpolomegaprime_new(ia,iz);
    lb = sum(vgridomega<nexta);
    lb(lb<=1) = 1;
    lb(lb>=pnumgridomega) = pnumgridomega-1;
    ub = lb+1;
    weightlb = (vgridomega(ub) - nexta)/(vgridomega(ub)-vgridomega(lb));
    weightlb(weightlb<0) = 0;
    weightlb(weightlb>1) = 1;
    weightub = 1-weightlb;

    mass = currentdist(ia,iz);
    for futureiz = 1:pnumgridz

        nextdist(lb,futureiz) = ...
            nextdist(lb,futureiz)...
            +mass*mtransz(iz,futureiz)*weightlb;

        nextdist(ub,futureiz) = ...
            nextdist(ub,futureiz)...
            +mass*mtransz(iz,futureiz)*weightub;

    end

end
end

%simulation routine;
error = max(abs(nextdist-currentdist),[],"all");
currentdist = nextdist;

end

%=========================
elseif eigvectormethod == 2
%=========================
%option2: eigenvector method (non-stochastic)
%=========================
%eigenvector method
mpolicy = zeros(pnumgridomega*pnumgridz,pnumgridomega*pnumgridz);

vlocationcombineda = kron(1:pnumgridomega,ones(1,pnumgridz));
vlocationcombinedz = kron(ones(size(vgridomega)),1:pnumgridz);

for iLocation = 1:pnumgridz*pnumgridomega

    ia = vlocationcombineda(iLocation);
    iz = vlocationcombinedz(iLocation);

    a = vgridomega(ia);
    nexta = mpolomegaprime_new(ia,iz);
    lb = sum(vgridomega<nexta);
    lb(lb<=0) = 1;
    lb(lb>=pnumgridomega) = pnumgridomega-1;
    ub = lb+1;
    weightlb = (vgridomega(ub) - nexta)/(vgridomega(ub)-vgridomega(lb));
    weightlb(weightlb<0) = 0;
    weightlb(weightlb>1) = 1;
    weightub = 1-weightlb;

    for izprime = 1:pnumgridz

        mpolicy(iLocation,:) = mpolicy(iLocation,:)+(vlocationcombineda==lb).*(vlocationcombinedz==izprime) * weightlb * mtransz(iz,izprime);
        mpolicy(iLocation,:) = mpolicy(iLocation,:)+(vlocationcombineda==ub).*(vlocationcombinedz==izprime) * weightub * mtransz(iz,izprime);

    end

end

mpolicy = sparse(mpolicy);
[currentDist0,~] = eigs(mpolicy',1);%eig(mPolicy');
currentDist0 = currentDist0(:,1)/sum(currentDist0(:,1));
currentdist = zeros(pnumgridomega,pnumgridz);

for iLocation = 1:pnumgridz*pnumgridomega

    ia = vlocationcombineda(iLocation);
    iz = vlocationcombinedz(iLocation);

    currentdist(ia,iz) = currentDist0(vlocationcombineda==ia & vlocationcombinedz==iz);

end
currentdist(currentdist<0) = 0;

end

%=========================  
% compute the equilibriun allocations
%=========================  
Omega           = sum(currentdist.*mgridomega,'all');
endoK           = sum(currentdist.*mpolkprime,'all');
endoB           = sum(currentdist.*mpolbprime,'all');
endoSupplyL     = sum(currentdist.*mpoln.*mgridz,'all');
Lambda          = sum(mlambda.*currentdist,'all'); % average lagrange multiplier
endoLambda      = sum(mlambda_new.*currentdist,'all'); % average lagrange multiplier
Phi             = sum(mphi.*currentdist,'all'); % average lagrange multiplier
endoPhi         = sum(mphi_new.*currentdist,'all'); % average lagrange multiplier
qnew            = sum(currentdist.*(mgridomega+mpoln.*mgridz*w-mpolc-mpolkprime+q*pBbar),'all')/pBbar;

%=========================  
% check the convergence and update the price
%=========================  
error2 = mean(abs(...
    [endoK - K;...
    endoSupplyL - supplyL;...
    qnew - q;...
    % Phi - endoPhi;...
    % endoB - 0;...
    ].^2), ...
    'all');

K           = weightold1.*K             + (1-weightold1).*endoK;
supplyL     = weightold2.*supplyL       + (1-weightold2).*endoSupplyL;
q           = weightold3.*q             + (1-weightold3).*qnew;

% report only spasmodically
if verbose == true && (floor((pnumiter_ge-1)/100) == (pnumiter_ge-1)/100) || error2<= tol_ge
%=========================  
% interim report
%=========================  

fprintf(' \n');
fprintf('market clearing results \n');
fprintf('max error: %.15f \n', error2);
fprintf('capital rent: %.15f \n', r);
fprintf('wage: %.15f \n', w);
fprintf('aggregate capital: %.15f \n', K);
fprintf('aggregate labor: %.15f \n', supplyL);
fprintf('bond price: %.15f \n', q);
fprintf('bond excess demand: %.15f \n', endoB);
% fprintf('bond excess demand: %.15f \n', endoB-pBbar);

% plot
close all;
figure;
plot(vgridomega,currentdist,'LineWidth',1.5);
xlim([min(vgridomega),max(vgridomega)]);
% legend('z_1','z_2','z_3','z_4','z_5','z_6','z_7',"FontSize",20);
saveas(gcf,'../figures/dist_ss.jpg');
% figure;
% plot(vgridomega,mpolomegaprime (:,1)); hold on;
% plot(vgridomega,mpolomegaprime (:,end));
% legend('Lowest z','Highest z','location','southeast',"FontSize",20);
% saveas(gcf,'../figures/policy.jpg');
% hold off;

pause(0.01);
toc;

%=========================
% save (mid)
%=========================
dir = '../solutions/WIP_ks1997endolabor_ss.mat';
save(dir);

end

pnumiter_ge = pnumiter_ge+1;

end

%=========================
%save
%=========================
dir = '../solutions/ks1997endolabor_ss.mat';
save(dir);
