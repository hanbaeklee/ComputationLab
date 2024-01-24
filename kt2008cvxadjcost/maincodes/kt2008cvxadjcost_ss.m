%% Solving the model in Khan and Thomas (2008)
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
eigvectormethod = false;        % eigenvector method on/off for stationary dist computation
verbose1        = false;        % vfi interim reports on/off
verbose2        = true;         % ge loop interim reports on/off

%=========================
% parameters - firms
%=========================
pnu             = 0.045; 
palpha          = 0.280;
pgamma          = 0.640;
pdelta          = 0.090;
pfc             = 0.420;
pcvx            = 0.890;

%=========================
% parameters - households
%=========================
pbeta           = 0.977;
peta            = 2.490;

%=========================
% numerical parameters - grids
%=========================
% finer grid near smaller capital (Maliar, Maliar, and Valli, 2010)
pnumgridk       = 50;
pgridkmin       = 1e-5;
pgridkmax       = 50;
x               = linspace(0,0.5,pnumgridk); 
y               = x.^5/max(x.^5);                   
vgridk          = pgridkmin+(pgridkmax-pgridkmin)*y;  
vgridkorig      = vgridk;

%=========================
% numerical parameters
%=========================
tol_ge          = 1e-10;
tol_vfi         = 1e-10;
tol_hist        = 1e-10;
weightold       = 0.95;
accinterval     = 30;
accstarting     = 30;

%=========================
% parameters - idiosyncratic productivity
%=========================
pmu = 0.0;
prho = 0.750;
psigma = 0.130;
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
vp = peta*ones(1);
mj = 0.01*sqrt(repmat(linspace(0.1,0.5,pnumgridk)',1,pnumgridz));
currentdist = ones(pnumgridk,pnumgridz)/(pnumgridk*pnumgridz);

%=========================    
% declare equilibrium objects
%=========================
mj0 = zeros(pnumgridk,pnumgridz);
mj1 = zeros(pnumgridk,pnumgridz);
futureval = zeros(pnumgridk,pnumgridz);
poleffdvd = zeros(pnumgridk,pnumgridz);
mkprimeopt = zeros(pnumgridk,pnumgridz);
poli = zeros(pnumgridk,pnumgridz);
polic = zeros(pnumgridk,pnumgridz);
poll = zeros(pnumgridk,pnumgridz);
poly = zeros(pnumgridk,pnumgridz);
polk = zeros(pnumgridk,pnumgridz);
polkc = zeros(pnumgridk,pnumgridz);
polthreshold = zeros(pnumgridk,pnumgridz);
polprob = zeros(pnumgridk,pnumgridz);

%%
%=========================    
% resume from the last one
%=========================    
% use the following line if you wish to start from where you stopped
% before.
% load '../solutions/WIP_kt2008_ss.mat';

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
    mj0 = mj*mtransz';
    z = mgridz;
    k = mgridk;                 
    l = ((pgamma)*z.*k.^(palpha)/w).^(1/(1-pgamma));
    
    % unconstrained optimum
    for iz = 1:pnumgridz
        
        kprimeopt = fnoptinvcvx(pbeta,pdelta,pcvx,vgridk,p,mj0(:,iz)',vgridk);
        low = sum(mgridkk <= kprimeopt',2);
        low(low<=0) = 1;
        low(low>=pnumgridk) = pnumgridk-1;
        high = low + 1;

        indxlow     = sub2ind(size(mgridkk),(1:pnumgridk)',low);
        indxhigh    = sub2ind(size(mgridkk),(1:pnumgridk)',high);
        
        weightdown  = (mgridkk(indxhigh) - kprimeopt')...
                    ./(mgridkk(indxhigh)-mgridkk(indxlow));
        weightup    = 1-weightdown;

        tempfutureval = mj0(low,iz).*weightdown +...
                        mj0(high,iz).*weightup;

        inv = (kprimeopt'-(1-pdelta)*vgridk');
        adjcost = pcvx/2*((inv./vgridk').^2).*vgridk';

        futureval(:,iz) =  ...
            - kprimeopt'*p - adjcost*p...
            + pbeta*tempfutureval;

        mkprimeopt(:,iz) = kprimeopt;

    end    
    
    % temporal profit
    tempprofit      = z.*k.^(palpha).*l.^(pgamma) ...
                    + (1-pdelta).*k - w*l;
    
    %-------------non-acceleration-------------%            
    if (floor((pnumiter_vfi-1)/accinterval) == (pnumiter_vfi-1)/accinterval || pnumiter_vfi<=accstarting)              

    % constrained optimum
    kprimeconst     = (mkprimeopt>(k*(1-pdelta+pnu))).*k*(1-pdelta+pnu) ...
                    + (mkprimeopt<(k*(1-pdelta-pnu))).*k*(1-pdelta-pnu) ...
                    + (mkprimeopt>=(k*(1-pdelta-pnu))).*(mkprimeopt<=(k*(1-pdelta+pnu))).*mkprimeopt;
    invconst        = kprimeconst-(1-pdelta)*k;
    adjcostconst    = pcvx/2*((invconst./k).^2).*k;
    
    tempfutureval   = interpn(mgridk, mgridz, mj0, ...
                                  kprimeconst, mgridz,'spline');    
    futurevalconst  = ...
                    - kprimeconst*p - adjcostconst*p...
                    + pbeta*tempfutureval ;
    
    % threshold
    tempxi          = (futureval - futurevalconst)*(1/w)*(1/p);
    xistar          = tempxi;
    xistar(xistar>pfc) = pfc;
    xistar(xistar<0) = 0;
    
    % optimizer
    optval          = tempprofit*p ...
                    + ((pfc-xistar)/pfc).*futurevalconst...
                    + (xistar/pfc).*(futureval-p*w*xistar/2);
    
    % adjustment cost
    inv             = mkprimeopt-(1-pdelta)*k;
    adjcost         = pcvx/2*((inv./k).^2).*k;
    
    % update
    mj1             = optval;
    poleffdvd       = z.*k.^(palpha).*l.^(pgamma) ...
                    - (inv + adjcost).*(xistar/pfc)...
                    - (invconst+adjcostconst).*((pfc-xistar)/pfc);
    poli            = mkprimeopt-(1-pdelta)*k;
    polic           = kprimeconst-(1-pdelta)*k;
    poll            = l+(xistar/pfc).*xistar/2;
    poly            = z.*k.^(palpha).*l.^(pgamma);
    polk            = mkprimeopt;
    polkc           = kprimeconst;                
    polthreshold    = xistar;
    polprob         = xistar/pfc;
    
    % measure the error (error is measured only in
    % the non-acceleration phase)
    error = max(abs(mj-mj1),[],'all');
    
    %-------------acceleration-------------%
    else
    
    kprimeopt       = polk;
    kprimeconst     = polkc;
    xistar          = polthreshold;
    invconst        = kprimeconst-(1-pdelta)*k;
    adjcostconst    = pcvx/2*((invconst./k).^2).*k;
    
    tempfutureval   = interpn(mgridk, mgridz, mj0, ...
                                  kprimeconst, mgridz,'spline');    
    futurevalconst  = ...
                    - kprimeconst*p - adjcostconst*p...
                    + pbeta*tempfutureval ;
    
    mj1             = tempprofit*p ...
                    + ((pfc-xistar)/pfc).*futurevalconst...
                    + (xistar/pfc).*(futureval-p*w*xistar/2);
    
    end                                                        
    
    % update the value function
    mj = mj1;
    
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
if eigvectormethod == true % eigenvector method

mpolicy = zeros(pnumgridk*pnumgridz,pnumgridk*pnumgridz);
vlocationcombinedk = kron(1:pnumgridk,ones(1,pnumgridz));
vlocationcombinedz = kron(ones(size(vgridk)),1:pnumgridz);

for ilocation = 1:pnumgridz*pnumgridk

ik = vlocationcombinedk(ilocation);
iz = vlocationcombinedz(ilocation);
xistar = polthreshold(ik,iz);

% lumpy
nextk = polk(ik,iz);
lb = sum(vgridk<nextk);
lb(lb<=0) = 1;
lb(lb>=pnumgridk) = pnumgridk-1;
ub = lb+1;
weightlb = (vgridk(ub) - nextk)/(vgridk(ub)-vgridk(lb));
weightlb(weightlb<0) = 0;
weightlb(weightlb>1) = 1;                    
weightub = 1-weightlb;

for izprime = 1:pnumgridz

    mpolicy(ilocation,:) = mpolicy(ilocation,:)+(vlocationcombinedk==lb).*(vlocationcombinedz==izprime) * weightlb * mtransz(iz,izprime) * (xistar/pfc);
    mpolicy(ilocation,:) = mpolicy(ilocation,:)+(vlocationcombinedk==ub).*(vlocationcombinedz==izprime) * weightub * mtransz(iz,izprime) * (xistar/pfc);

end

% non-lumpy
nextk = polkc(ik,iz);
lb = sum(vgridk<nextk);
lb(lb<=0) = 1;
lb(lb>=pnumgridk) = pnumgridk-1;
ub = lb+1;
weightlb = (vgridk(ub) - nextk)/(vgridk(ub)-vgridk(lb));
weightlb(weightlb<0) = 0;
weightlb(weightlb>1) = 1;                    
weightub = 1-weightlb;

for izprime = 1:pnumgridz

    mpolicy(ilocation,:) = mpolicy(ilocation,:)+(vlocationcombinedk==lb).*(vlocationcombinedz==izprime) * weightlb * mtransz(iz,izprime) * ((pfc-xistar)/pfc);
    mpolicy(ilocation,:) = mpolicy(ilocation,:)+(vlocationcombinedk==ub).*(vlocationcombinedz==izprime) * weightub * mtransz(iz,izprime) * ((pfc-xistar)/pfc);

end

end
mpolicy_trans = sparse(mpolicy');

[currentdist0,~] = eigs(mpolicy_trans,1);%eig(mPolicy');
currentdist0 = currentdist0(:,1)/sum(currentdist0(:,1));
currentdist = zeros(pnumgridk,pnumgridz);

for ilocation = 1:pnumgridz*pnumgridk

    ik = vlocationcombinedk(ilocation);
    iz = vlocationcombinedz(ilocation);

    currentdist(ik,iz) = currentdist0(vlocationcombinedk==ik & vlocationcombinedz==iz);

end

else % histogram-evolving method    

error = 10;
while error>tol_hist

% empty distribution to save the next distribution  
nextdist = zeros(size(currentdist));

for iz = 1:pnumgridz

% policies   
nextkall = polk(:,iz); 
nextkcall = polkc(:,iz);
xistar = polthreshold(:,iz);
   
% lumpy
nextk = nextkall;
lb = sum(mgridkk<nextk,2);
lb(lb<=0) = 1;
lb(lb>=pnumgridk) = pnumgridk-1;
ub = lb+1;
weightlb = (vgridk(ub) - nextk')./(vgridk(ub)-vgridk(lb));
weightlb(weightlb<0) = 0;
weightlb(weightlb>1) = 1;                    
weightub = 1-weightlb;

mass = currentdist(:,iz) .* (xistar/pfc);

for futureiz = 1:pnumgridz
for ik = 1:pnumgridk    
    nextdist(lb(ik),futureiz) = ...
        nextdist(lb(ik),futureiz)...
        +mass(ik).*mtransz(iz,futureiz).*weightlb(ik);

    nextdist(ub(ik),futureiz) = ...
        nextdist(ub(ik),futureiz)...
        +mass(ik).*mtransz(iz,futureiz).*weightub(ik);
end
end 

% nonlumpy
nextk = nextkcall;
lb = sum(mgridkk<nextk,2);
lb(lb<=0) = 1;
lb(lb>=pnumgridk) = pnumgridk-1;
ub = lb+1;
weightlb = (vgridk(ub) - nextk')./(vgridk(ub)-vgridk(lb));
weightlb(weightlb<0) = 0;
weightlb(weightlb>1) = 1;                    
weightub = 1-weightlb;

mass = currentdist(:,iz) .* ((pfc-xistar)/pfc);

for futureiz = 1:pnumgridz
for ik = 1:pnumgridk    
    nextdist(lb(ik),futureiz) = ...
        nextdist(lb(ik),futureiz)...
        +mass(ik).*mtransz(iz,futureiz).*weightlb(ik);

    nextdist(ub(ik),futureiz) = ...
        nextdist(ub(ik),futureiz)...
        +mass(ik).*mtransz(iz,futureiz).*weightub(ik);
end
end 
           
end

error = max(abs(nextdist-currentdist),[],'all');
currentdist = nextdist;

end
end

%=========================  
% 4. compute the equilibrium allocations
%=========================  
% fraction of firms with (i/k > 0.2)
lumpy = ((poli./repmat(vgridk',1,pnumgridz))>0.2).*polprob ...
       +((polic./repmat(vgridk',1,pnumgridz))>0.2).*(1-polprob);

% fraction of firms with (i > 0.0)
lumpy0 = ((poli./repmat(vgridk',1,pnumgridz))>0.0).*polprob ...
       +((polic./repmat(vgridk',1,pnumgridz))>0.0).*(1-polprob);

% fraction of firms with (i/k < -0.2)
neglumpy = ((poli./repmat(vgridk',1,pnumgridz))<-0.2).*polprob ...
          +((polic./repmat(vgridk',1,pnumgridz))<-0.2).*(1-polprob);

% fraction of firms with (|i/k| < 0.01)
inaction = (abs(poli./repmat(vgridk',1,pnumgridz))<0.01).*(polprob)...
          +(abs(polic./repmat(vgridk',1,pnumgridz))<0.01).*(1-polprob);

% fraction of firms with (i/k > 0.01)
posinv = ((poli./repmat(vgridk',1,pnumgridz))>0.01).*(polprob)...
        +((polic./repmat(vgridk',1,pnumgridz))>0.01).*(1-polprob);

% fraction of firms with (i/k < -0.01)
neginv = ((poli./repmat(vgridk',1,pnumgridz))<-0.01).*(polprob)...
        +((polic./repmat(vgridk',1,pnumgridz))<-0.01).*(1-polprob);

eq.p = vp;
eq.w = vw;  

eq.k = sum(repmat(vgridk',1,pnumgridz).*currentdist,'all');
eq.i = sum((poli.*(polprob)+polic.*(1-polprob)).*currentdist,'all');
eq.ik = sum(((poli.*(polprob)+polic.*(1-polprob))./repmat(vgridk',1,pnumgridz)).*currentdist,'all');
eq.lumpy = sum(lumpy.*currentdist,'all');
eq.ilumpy = sum(lumpy0.*(poli.*(polprob)+polic.*(1-polprob)).*currentdist,'all')/sum(lumpy0.*currentdist,'all');
eq.l = sum(poll.*currentdist,'all');
eq.y = sum(poly.*currentdist,'all');
eq.c = sum(poleffdvd.*currentdist,'all');

%=========================  
% general equilibrium routine
%=========================  
eq.p_implied = 1./eq.c;
errorp = eq.p_implied-eq.p;        

%=========================      
% 5. simulated method of moments
%=========================  
moments.kt.inaction = sum(inaction.*currentdist,'all');
moments.kt.posspike = sum(lumpy.*currentdist,'all');
moments.kt.negspike = sum(neglumpy.*currentdist,'all');
moments.kt.posinv = sum(posinv.*currentdist,'all');
moments.kt.neginv = sum(neginv.*currentdist,'all');

moments.wb.ikratio = sum(((poli.*(polprob)+polic.*(1-polprob))./repmat(vgridk',1,pnumgridz)).*currentdist,'all');
moments.wb.ikratio_sd = sqrt(sum((((poli./repmat(vgridk',1,pnumgridz))-moments.wb.ikratio).^2).*currentdist.*polprob,'all')...
                  +sum((((polic./repmat(vgridk',1,pnumgridz))-moments.wb.ikratio).^2).*currentdist.*(1-polprob),'all'));
moments.wb.spikeratio = eq.lumpy;
moments.wb.posinvrate = 1-eq.lumpy;

moments.mktclr = errorp;

% for potential calibration/estimation
vtarget = [moments.wb.spikeratio;moments.wb.ikratio;moments.wb.ikratio_sd;moments.mktclr];

%=========================  
% 6. update the price
%=========================  
oldGridp = vp; % to save the prior price vector
vp = vp.*weightold+eq.p_implied.*(1-weightold);
error2 = mean(abs(errorp).^2);

% report only spasmodically
if verbose2 == true && (floor((pnumiter_ge-1)/1) == (pnumiter_ge-1)/1) || error2<= tol_ge
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
fprintf('Khan and Thomas (2008) moments\n');  
disp(moments.kt);  
fprintf('Winberry (2021) moments\n');  
disp(moments.wb);  
            
toc;
%=========================  
% save (mid)
%=========================  
save '../solutions/WIP_kt2008cvxadjcost_ss.mat';    
end

pnumiter_ge = pnumiter_ge+1;

end

%%
%=========================  
% save (final)
%=========================  
save '../solutions/kt2008cvxadjcost_ss.mat';

