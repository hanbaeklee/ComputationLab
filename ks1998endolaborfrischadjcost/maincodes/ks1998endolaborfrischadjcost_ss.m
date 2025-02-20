%% Krusell and Smith (1998) with endogenous labor supply and convex adjustment cost
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
eigvectormethod1 = 2;            % eigenvector method on/off for labor supply computation
eigvectormethod2 = 2;            % eigenvector method on/off for stationary dist computation
verbose          = true;         % ge loop interim reports on/off

%=========================
% parameters
%=========================
palpha      = 0.33; 
pbeta       = 0.99; 
pdelta      = 0.025;
prho        = 0.90;
psigma      = 0.05;
pfrisch     = 1.00;
peta        = 7.60;
pmu         = 0.60;

%=========================
% numerical parameters - grids
%=========================
%idiosyncratic income shock
pnumgridz   = 7;
[mtransz, vgridz] = ...
fnTauchen(prho, 0, psigma^2, pnumgridz, 3);
vgridz = exp(vgridz');            
if pnumgridz ==1
    mtransz = 1;
end

%wealth grid
pnumgrida   = 100;

%finer grid near smaller wealth (Maliar, Maliar, and Valli, 2010)
vgridamin   = 1.0;
vgridamax   = 200;
x           = linspace(0,0.5,pnumgrida);
y           = x.^7/max(x.^7);
vgrida      = vgridamin+(vgridamax-vgridamin)*y;

%=========================
% numerical parameters
%=========================
tol_labor       = 1e-10;
tol_ge          = 1e-8;
tol_pfi         = 1e-10;
tol_hist        = 1e-10;
weightold1      = 0.9000;
weightold2      = 0.9000;
weightold3      = 0.9000;
weightold4      = 0.9000;

%=========================    
% extra setup
%=========================    
% grids for matrix operations
mgrida      = repmat(vgrida',1,pnumgridz);
mgridz      = repmat(vgridz,pnumgrida,1);
mgridaa     = repmat(vgrida,pnumgrida,1);

%=========================    
% initial guess
%=========================    
K           = 7.50;
supplyL     = 0.33;
currentdist = ones(pnumgrida,pnumgridz)/(pnumgrida*pnumgridz);

%=========================    
% declare equilibrium objects
%=========================
mpolc       = repmat(0.01*vgrida',1,pnumgridz);
mpolaprime  = ones(size(mpolc));
mpolaprime_new = zeros(size(mpolc));
mlambda     = zeros(size(mpolc));
mlambda_new = zeros(size(mpolc));


%%
%=========================    
% resume from the last one
%=========================    
% use the following line if you wish to start from where you stopped
% before.
load '../solutions/WIP_ks1998endolaborfrischadjcost_ss.mat';

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

mexpectation = 0;
for izprime = 1:pnumgridz
    
    zprime = vgridz(izprime);
    mpolaprimeprime = interp1(vgrida',squeeze(mpolaprime(:,izprime)),...
                    squeeze(mpolaprime),"linear","extrap"); 
    rprime = r;
    wprime = w;
    psi2 = -(pmu/2)*((mpolaprimeprime./mpolaprime).^2-1);

    mprime = ((1+rprime)*mpolaprime - (pmu/2).*((mpolaprimeprime-mpolaprime)./mpolaprime).^2.*mpolaprime - mpolaprimeprime)/(wprime*zprime); 
    nprime = (-peta*mprime + sqrt((peta*mprime).^2+4*peta))/(2*peta);
    cprime = wprime.*zprime*nprime + (1+rprime).*squeeze(mpolaprime) - mpolaprimeprime  ...
           - (pmu/2).*((mpolaprimeprime-mpolaprime)./mpolaprime).^2.*mpolaprime;
    cprime(cprime<=0) = 1e-10 ;
    muprime = 1./cprime;                
    mexpectation =  mexpectation + repmat(mtransz(:,izprime)',pnumgrida,1).*(1+rprime - psi2).*muprime;
    
end

mexpectation = pbeta*mexpectation;
c = (1+pmu*(mpolaprime-mgrida)./mgrida)./(mexpectation+mlambda);
mpoln = (w.*mgridz./(peta*c)).^pfrisch;
mlambda_new = (1+pmu*(mpolaprime-mgrida)./mgrida)./(w.*mgridz.*mpoln + (1+r).*mgrida - (pmu/2).*((mpolaprime-mgrida)./mgrida).^2.*mgrida ...
            - mpolaprime) - mexpectation;
mpolaprime_new = w.*mgridz.*mpoln + (1+r).*vgrida' - c - (pmu/2).*((mpolaprime-mgrida)./mgrida).^2.*mgrida;

mlambda_new(mpolaprime_new>vgridamin) = 0;
mpolaprime_new(mpolaprime_new<=vgridamin) = vgridamin;
mpolc = c;

% iteration routine
error = max(abs(mpolaprime-mpolaprime_new),[],"all");


%%
%=========================
if eigvectormethod2 == 1
%=========================
%option1: histogram-evolving method (non-stochastic)
%=========================
error = 10;
while error>tol_hist
% empty distribution to save the next distribution  
nextdist = zeros(size(currentdist));

for iz = 1:pnumgridz
for ia = 1:pnumgrida

    a = vgrida(ia);
    nexta = mpolaprime_new(ia,iz);
    lb = sum(vgrida<nexta);
    lb(lb<=1) = 1;
    lb(lb>=pnumgrida) = pnumgrida-1;
    ub = lb+1;
    weightlb = (vgrida(ub) - nexta)/(vgrida(ub)-vgrida(lb));
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
elseif eigvectormethod2 == 2
%=========================
%option2: eigenvector method (non-stochastic)
%=========================
%eigenvector method
mpolicy = zeros(pnumgrida*pnumgridz,pnumgrida*pnumgridz);

vlocationcombineda = kron(1:pnumgrida,ones(1,pnumgridz));
vlocationcombinedz = kron(ones(size(vgrida)),1:pnumgridz);

for iLocation = 1:pnumgridz*pnumgrida

    ia = vlocationcombineda(iLocation);
    iz = vlocationcombinedz(iLocation);

    a = vgrida(ia);
    nexta = mpolaprime_new(ia,iz);
    lb = sum(vgrida<nexta);
    lb(lb<=0) = 1;
    lb(lb>=pnumgrida) = pnumgrida-1;
    ub = lb+1;
    weightlb = (vgrida(ub) - nexta)/(vgrida(ub)-vgrida(lb));
    weightlb(weightlb<0) = 0;
    weightlb(weightlb>1) = 1;
    weightub = 1-weightlb;

    for izprime = 1:pnumgridz

        mpolicy(iLocation,:) = mpolicy(iLocation,:)+(vlocationcombineda==lb).*(vlocationcombinedz==izprime) * weightlb * mtransz(iz,izprime);
        mpolicy(iLocation,:) = mpolicy(iLocation,:)+(vlocationcombineda==ub).*(vlocationcombinedz==izprime) * weightub * mtransz(iz,izprime);

    end

end

[currentDist0,~] = eigs(mpolicy',1);%eig(mPolicy');
currentDist0 = currentDist0(:,1)/sum(currentDist0(:,1));
currentdist = zeros(pnumgrida,pnumgridz);

for iLocation = 1:pnumgridz*pnumgrida

    ia = vlocationcombineda(iLocation);
    iz = vlocationcombinedz(iLocation);

    currentdist(ia,iz) = currentDist0(vlocationcombineda==ia & vlocationcombinedz==iz);

end
currentdist(currentdist<0) = 0;

end

%=========================  
% compute the equilibriun allocations
%=========================  
marginalDista   = sum(currentdist,2);
endoK           = vgrida*marginalDista;
Lambda          = sum(mlambda.*currentdist,'all'); % average lagrange multiplier
endoLambda      = sum(mlambda_new.*currentdist,'all'); % average lagrange multiplier
endoSupplyL     = sum(currentdist.*mpoln.*mgridz,'all');

%=========================  
% check the convergence and update the price
%=========================  
error2 = mean(abs(...
    [endoK - K;...
    endoSupplyL - supplyL;...
    Lambda - endoLambda]), ...
    'all');

K           = weightold1.*K             + (1-weightold1).*endoK;
supplyL     = weightold2.*supplyL       + (1-weightold2).*endoSupplyL;
mlambda     = weightold3.*mlambda       + (1-weightold3).*mlambda_new;
mpolaprime  = weightold4.*mpolaprime    + (1-weightold4).*mpolaprime_new;

% report only spasmodically
if verbose == true && (floor((pnumiter_ge-1)/200) == (pnumiter_ge-1)/200) || error2<= tol_ge
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

% plot
close all;
figure;
plot(vgrida,currentdist,'LineWidth',1.5);
legend('z_1','z_2','z_3','z_4','z_5','z_6','z_7',"FontSize",20);
saveas(gcf,'../figures/dist_ss.jpg');
figure;
plot(vgrida,mpolaprime (:,1)); hold on;
plot(vgrida,mpolaprime (:,end));
legend('Lowest z','Highest z','location','southeast',"FontSize",20);
saveas(gcf,'../figures/policy.jpg');
hold off;

pause(0.01);
toc;

%=========================
% save (mid)
%=========================
dir = '../solutions/WIP_ks1998endolaborfrischadjcost_ss.mat';
save(dir);

end

pnumiter_ge = pnumiter_ge+1;

end

%=========================
%save
%=========================
dir = '../solutions/ks1998endolaborfrischadjcost_ss.mat';
save(dir);
