%% Solving the model in Krusell and Smith (1998)
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
palpha     = 0.36; 
pbeta      = 0.99; 
pdelta     = 0.025;
ptheta     = 1/2.9;

%=========================
% numerical parameters - grids
%=========================
%idiosyncratic income shock
pnumgridz   = 2;
% vgridz      = [0.00,1.00];
vgridz      = [0,1.00];
pnumgridA   = 2;
mtransz     = zeros(pnumgridz,pnumgridz);
mtransz0    = [0.525,0.350,0.03125,0.09375;...
               0.035,0.84,0.0025,0.1225;...
               0.09375,0.03125,0.292,0.583;...
               0.0099,0.1151,0.0245,0.8505];    
mtransA = [sum(sum(mtransz0(1:2,1:2))),sum(sum(mtransz0(1:2,3:4))); ...
           sum(sum(mtransz0(3:4,1:2))),sum(sum(mtransz0(3:4,3:4)))];
mtransA = mtransA/(2*pnumgridA);
mtransz0(1,1:2) = mtransz0(1,1:2)/sum(mtransz0(1,1:2),'all');
mtransz0(2,1:2) = mtransz0(2,1:2)/sum(mtransz0(2,1:2),'all');
mtransz0(1,3:4) = mtransz0(1,3:4)/sum(mtransz0(1,3:4),'all');
mtransz0(2,3:4) = mtransz0(2,3:4)/sum(mtransz0(2,3:4),'all');
mtransz0(3,1:2) = mtransz0(3,1:2)/sum(mtransz0(3,1:2),'all');
mtransz0(4,1:2) = mtransz0(4,1:2)/sum(mtransz0(4,1:2),'all');
mtransz0(3,3:4) = mtransz0(3,3:4)/sum(mtransz0(3,3:4),'all');
mtransz0(4,3:4) = mtransz0(4,3:4)/sum(mtransz0(4,3:4),'all');

% krusell and smith (1998) assumes the idiosyncratic shock transition
% probability depends upon the aggregate shock realization. Thus, the 
% stationary equilibrium calculation needs an extra step to compute the
% stationary transition probability as follows:

%uu
mtransz(1,1)= mtransz0(1,1)*mtransA(1,1)...
             +mtransz0(1,3)*mtransA(1,2)...
             +mtransz0(3,1)*mtransA(2,1)...
             +mtransz0(3,3)*mtransA(2,2);
%ue
mtransz(1,2)= mtransz0(1,2)*mtransA(1,1)...
             +mtransz0(1,4)*mtransA(1,2)...
             +mtransz0(3,2)*mtransA(2,1)...
             +mtransz0(3,4)*mtransA(2,2);
%eu
mtransz(2,1)= mtransz0(2,1)*mtransA(1,1)...
             +mtransz0(2,3)*mtransA(1,2)...
             +mtransz0(4,1)*mtransA(2,1)...
             +mtransz0(4,3)*mtransA(2,2);
%ee
mtransz(2,2)= mtransz0(2,2)*mtransA(1,1)...
             +mtransz0(2,4)*mtransA(1,2)...
             +mtransz0(4,2)*mtransA(2,1)...
             +mtransz0(4,4)*mtransA(2,2);

%wealth grid
pnumgrida   = 100;

%finer grid near smaller wealth (Maliar, Maliar, and Valli, 2010)
vgridamin   = 1e-20;%0;
vgridamax   = 100;
x           = linspace(0,0.5,pnumgrida);
y           = x.^7/max(x.^7);
vgrida      = vgridamin+(vgridamax-vgridamin)*y;

%=========================
% numerical parameters
%=========================
tol_labor       = 1e-10;
tol_ge          = 1e-10;
tol_pfi         = 1e-10;
tol_hist        = 1e-10;
weightold1      = 0.9990;
weightold2      = 0.9990;
weightold3      = 0.8000;
weightold4      = 0.8000;
weightold5      = 0.8000;

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
K           = 11;
supplyL     = 0.28;%0.33;
currentdist = ones(pnumgrida,pnumgridz)/(pnumgrida*pnumgridz);

%=========================    
% declare equilibrium objects
%=========================
mpolc       = repmat(0.01*vgrida',1,pnumgridz);
mpoln       = ones(size(mpolc))*supplyL;
mpolaprime  = zeros(size(mpolc));
mpolaprime_new = zeros(size(mpolc));
mlambda     = zeros(size(mpolc));
mlambda_new = zeros(size(mpolc));


%%
%=========================    
% resume from the last one
%=========================    
% use the following line if you wish to start from where you stopped
% before.
load '../solutions/WIP_ks1998endolabor_ss.mat';

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
pnumiter_vfi = 1;
while error>tol_pfi
mexpectation = 0;    
for ishockprime = 1:pnumgridz
    
    izprime = ishockprime;
    zprime = vgridz(izprime);
    
    mpolaprimeprime = interp1(vgrida',squeeze(mpolaprime(:,izprime)),...
        mpolaprime,"linear","extrap"); 
    mpolnprime = interp1(vgrida',squeeze(mpoln(:,izprime)),...
        mpolaprime,"linear","extrap"); 
    rprime = r;
    wprime = w;
    cprime = wprime.*zprime.*mpolnprime + (1+rprime).*mpolaprime - mpolaprimeprime;
    cprime(cprime<=1e-10) = 1e-10;
    muprime = 1./cprime;                
    mexpectation =  mexpectation + repmat(mtransz(:,izprime)',pnumgrida,1).*(1+rprime).*muprime;
    
end

c = 1./(pbeta*mexpectation + mlambda);
mlambda_new = 1./(w.*mgridz.*mpoln + (1+r).*mgrida - mpolaprime) - pbeta*mexpectation;
mlambda_new(mlambda_new<0) = 0;
mpoln_new = 1-(1-ptheta).*mpolc./(w.*mgridz.*ptheta);
mpoln_new(mpoln_new>1) = 1;
mpoln_new(mpoln_new<0) = 0;
mpolaprime_new = w.*mgridz.*mpoln + (1+r).*mgrida - c;

mlambda_new(mpolaprime_new>vgridamin) = 0;
mpolaprime_new(mpolaprime_new<=vgridamin) = vgridamin;
mpolc = c;
mpolc(mpolc<0) = 0;

% iteration routine
error = max(abs(mpolaprime-mpolaprime_new),[],"all");
pnumiter_vfi = pnumiter_vfi+1;

mlambda     = weightold3.*mlambda       + (1-weightold3).*mlambda_new;
mpolaprime  = weightold4.*mpolaprime    + (1-weightold4).*mpolaprime_new;
mpoln       = weightold5.*mpoln         + (1-weightold5).*mpoln_new;

end

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

[currentdist0,~] = eigs(mpolicy',1);%eig(mPolicy');
currentdist0 = currentdist0(:,1)/sum(currentdist0(:,1));
currentdist = zeros(pnumgrida,pnumgridz);

for iLocation = 1:pnumgridz*pnumgrida

    ia = vlocationcombineda(iLocation);
    iz = vlocationcombinedz(iLocation);

    currentdist(ia,iz) = currentdist0(vlocationcombineda==ia & vlocationcombinedz==iz);

end
currentdist(currentdist<0) = 0;

end

%=========================  
% compute the equilibriun allocations
%=========================  
marginaldista   = sum(currentdist,2);
endoK           = vgrida*marginaldista;
Lambda          = sum(mlambda.*currentdist,'all'); % average lagrange multiplier
endoLambda      = sum(mlambda_new.*currentdist,'all'); % average lagrange multiplier
endoSupplyL     = sum(currentdist.*mpoln.*mgridz,'all');

%=========================  
% check the convergence and update the price
%=========================  
error2 = mean(abs(...
    [endoK - K;...
    endoSupplyL - supplyL]), ...
    'all');

K           = weightold1.*K             + (1-weightold1).*endoK;
supplyL     = weightold2.*supplyL       + (1-weightold2).*endoSupplyL;

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
% legend('z_1','z_2','z_3','z_4','z_5','z_6','z_7',"FontSize",20);
legend('z_1','z_2',"FontSize",20);
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
dir = '../solutions/WIP_ks1998endolabor_ss.mat';
save(dir);

end

pnumiter_ge = pnumiter_ge+1;

end

%=========================
%save
%=========================
dir = '../solutions/ks1998endolabor_ss.mat';
save(dir);
