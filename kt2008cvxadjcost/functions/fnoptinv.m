function [kprimeopt] = fnOptInv(pbeta,pdelta,pgammaI,vgridkorig,p,matjsprime)

    plength = length(vgridkorig);    
    vgridkorigup = [vgridkorig(2:plength),vgridkorig(plength)]';
    vgridkorigdn = [vgridkorig(1),vgridkorig(1:(plength-1))]';
    gap = vgridkorigup-vgridkorigdn;
    
    futureval1 = pbeta*[matjsprime(2:plength);matjsprime(plength)];
    futureval0 = pbeta*[matjsprime(1);matjsprime(1:(plength-1))];
    
    futurevalderiv = (futureval1-futureval0)./(gap);
    vfoc = futurevalderiv - p*pgammaI;
    
    if length(unique(vfoc))~=length(vfoc)
        vfoc = vfoc-linspace(0,0.0000001,length(vfoc))';
    end
    
    % interpolation
    low = sum(vfoc>0);
    low(low<1) = 1;
    low(low>=(length(vfoc)-1)) = length(vfoc)-1;
    high = low + 1;
        
    weightLow = (vfoc(high) - 0)/(vfoc(high)-vfoc(low));
    kprimeopt = weightLow*vgridkorig(low)+(1-weightLow)*vgridkorig(high);    
    
    lb = min(vgridkorig);
    ub = max(vgridkorig);
    
    kprimeopt(kprimeopt<lb) = lb;
    kprimeopt(kprimeopt>ub) = ub;
    
end