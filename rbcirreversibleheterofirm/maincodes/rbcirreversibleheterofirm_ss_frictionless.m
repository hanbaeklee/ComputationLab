%% An RBC model with heterogeneous firms and irreversible investments
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
verbose1        = false;         % vfi interim reports on/off
verbose2        = true;          % ge loop interim reports on/off

%=========================
% parameters - firms
%=========================
palpha          = 0.256;
pgamma          = 0.640;
pdelta          = 0.069;
% pphi            = 0.975;
pphi            = -9999;

%=========================
% parameters - households
%=========================
pbeta           = 0.977;
peta            = 2.40;

%=========================
% numerical parameters - grids
%=========================
% finer grid near smaller capital (Maliar, Maliar, and Valli, 2010)
pnumgridk       = 100;
pgridkmin       = 1e-5;
pgridkmax       = 10;
x               = linspace(0,0.5,pnumgridk); 
y               = x.^3/max(x.^3);                   
vgridk          = pgridkmin+(pgridkmax-pgridkmin)*y;  
vgridkorig      = vgridk;

%=========================
% numerical parameters
%=========================
tol_ge          = 1e-10;
tol_vfi         = 1e-10;
tol_hist        = 1e-10;
weightold       = 0.90;
accinterval     = 30;
accstarting     = 10;

%=========================
% parameters - idiosyncratic productivity
%=========================
pmu = 0.0;
prho = 0.859;
psigma = 0.022;
pnumgridz = 11;
[mtransz, vgridz] = ...
fnTauchen(prho, pmu, psigma^2, pnumgridz, 3);
vgridz = exp(vgridz');                
if pnumgridz ==1
    mtransz = 1;
end

%=========================    
% extra setup
%=========================    
% grids for matrix operations
mgridk          = repmat(vgridk',1,pnumgridz);
mgridz          = repmat(vgridz,pnumgridk,1);
mgridkk         = repmat(vgridk,pnumgridk,1);

%=========================    
% initial guess
%=========================    
vp = 2.5*ones(1);
threshold = 0.10;
mjderiv = vp*(1-pgamma)*(pgamma/(peta/vp)).^(pgamma/(1-pgamma)).*mgridz.^(1/(1-pgamma))...
        .* (palpha/(1-pgamma)).*mgridk.^(palpha/(1-pgamma)-1) ...
        + vp*(1-pdelta);
currentdist = ones(pnumgridk,pnumgridz)/(pnumgridk*pnumgridz);

%=========================    
% declare equilibrium objects
%=========================
mj1deriv = zeros(pnumgridk,pnumgridz);
poleffdvd = zeros(pnumgridk,pnumgridz);
poli = zeros(pnumgridk,pnumgridz);
poll = zeros(pnumgridk,pnumgridz);
poly = zeros(pnumgridk,pnumgridz);
polk = zeros(pnumgridk,pnumgridz);
mlambda = zeros(pnumgridk,pnumgridz);

mkprimeopt = zeros(pnumgridk,pnumgridz);
futurevalderiv = zeros(pnumgridk,pnumgridz);

%%
%=========================    
% resume from the last one
%=========================    
% use the following line if you wish to start from where you stopped
% before.
% load '../solutions/WIP_rbcirreversibleheterofirm_ss_frictionless.mat';

%%
%=========================    
% outer loop
%=========================    
tic;
error2 = 10;
pnumiter_ge = 1;

while error2>tol_ge

% macro setup: wage
vw = peta./vp;    

% auxiliary
pnumiter_vfi = 1;
error = 10;    

%=========================    
% 2. vfi
%=========================  
while error>tol_vfi
     
    w = vw;
    p = vp;
    
    % expected future value
    mj0deriv = mjderiv*mtransz';
    z = mgridz;
    k = mgridk;                 
    l = ((pgamma)*z.*k.^(palpha)/w).^(1/(1-pgamma));
   
    for iz = 1:pnumgridz
        
        tempmat = mj0deriv(:,iz)';
        kprimeopt = fnoptinvirreversible(pbeta,vgridk,p,tempmat,mlambda(:,iz));
        kprimeopt = (threshold + (1-pdelta)*vgridk ).*(kprimeopt< threshold+(1-pdelta) *vgridk) ...
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

        tempfuturevalderiv = mj0deriv(low,iz).*weightdown + ...
                             mj0deriv(high,iz).*weightup;
        lambda =  p - pbeta*tempfuturevalderiv;
        lambda(lambda<0) = 0;
        
        mkprimeopt(:,iz) = kprimeopt;
        mlambda(:,iz) = lambda;
        %mlambda(:,iz) = 0;

    end    
    
    % updated value
    optvalderiv     = p*(1-pgamma)*(pgamma/w).^(pgamma/(1-pgamma)).*z.^(1/(1-pgamma))...
                    .* (palpha/(1-pgamma)).*k.^(palpha/(1-pgamma)-1) ...
                    + (1-pdelta).*(p-mlambda);
    % update
    mj1deriv        = optvalderiv;
    poleffdvd       = z.*k.^(palpha).*l.^(pgamma) - (mkprimeopt-(1-pdelta)*k);
    poli            = mkprimeopt-(1-pdelta)*k;
    poll            = l;
    poly            = z.*k.^(palpha).*l.^(pgamma);
    polk            = mkprimeopt;

    % measure the error (error is measured only in
    % the non-acceleration phase)
    error = max(abs(mjderiv-mj1deriv),[],'all');
                                                
    % update the value function
    mjderiv = mj1deriv;
    
    if verbose1 == true && (floor((pnumiter_vfi-1)/30) == (pnumiter_vfi-1)/30)
    dispthis=['iteration is in process : ', num2str(pnumiter_vfi),' X ',num2str(error)];
    disp(dispthis)
    toc;   
    end
    
    pnumiter_vfi = pnumiter_vfi+1;

end
    
%=========================    
% 3. non-stochastic simulation
%=========================       
% histogram-evolving method    
error = 10;
while error>tol_hist

% empty distribution to save the next distribution  
nextdist = zeros(size(currentdist));

for iz = 1:pnumgridz
for ik = 1:pnumgridk

nextk = polk(ik,iz);
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

error = max(abs(nextdist-currentdist),[],'all');
currentdist = nextdist;

end

%=========================  
% 4. compute the equilibrium allocations
%=========================  
eq.p = vp;
eq.w = vw;  

eq.k = sum(repmat(vgridk',1,pnumgridz).*currentdist,'all');
eq.i = sum(poli.*currentdist,'all');
eq.ik = sum(poli./repmat(vgridk',1,pnumgridz).*currentdist,'all');
eq.l = sum(poll.*currentdist,'all');
eq.y = sum(poly.*currentdist,'all');
eq.c = sum(poleffdvd.*currentdist,'all');
eq.threshold = eq.i*pphi;

%=========================  
% general equilibrium routine
%=========================  
eq.p_implied = 1./eq.c;
errorp = eq.p_implied-eq.p;        
errorthreshold = eq.threshold-threshold;        

%=========================  
% 5. update the price
%=========================  
oldGridp    = vp; % to save the prior price vector
vp          = vp.*weightold         + eq.p_implied.*(1-weightold);
threshold   = threshold*weightold   + eq.threshold.*(1-weightold);
error2      = mean([abs(errorp),abs(errorthreshold)].^2);

% report only spasmodically
if verbose2 == true && (floor((pnumiter_ge-1)/10) == (pnumiter_ge-1)/10) || error2<= tol_ge
%=========================  
% interim report
%=========================  
subplot(1,2,1)
for iz = 1:pnumgridz
    plot(vgridk,currentdist(:,iz)./sum(currentdist(:,iz)));hold on;
end
xlim([min(vgridk),max(vgridk)]);
hold off;
title("productivity-specific capital dist.")
subplot(1,2,2)       
plot(vgridk,sum(currentdist,2));hold on;
xlim([min(vgridk),max(vgridk)]);
hold off;
title("marginal capital dist.")
legend("Marginal capital distribution","location","northeast");
pause(0.01);
location = ['../figures/dist_ss.pdf'];
saveas(gcf, location);

fprintf(' \n');
fprintf('Market Clearing Results \n');
fprintf('Max error: %.10f \n', error2);  
fprintf('Price: %.10f \n', vp);  
fprintf('Labor: %.10f \n', eq.l); 
fprintf('Price: %.10f \n', vp);  
fprintf(' \n');  
            
toc;
%=========================  
% save (mid)
%=========================  
save '../solutions/WIP_rbcirreversibleheterofirm_ss_frictionless.mat';    
end

pnumiter_ge = pnumiter_ge+1;

end

%%
%=========================  
% save (final)
%=========================  
save '../solutions/rbcirreversibleheterofirm_ss_frictionless.mat';

