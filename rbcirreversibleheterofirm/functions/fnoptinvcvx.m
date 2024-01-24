function [kprimeopt] = fnoptinvcvx(pbeta,pdelta,pcvx,vgridk,p,mj,k)

    plength = length(vgridk);
    
    vgridkup = [vgridk(2:plength),vgridk(plength)];
    vgridkdn = [vgridk(1),vgridk(1:(plength-1))];
    gap = vgridkup-vgridkdn;
    
    futureval1 = pbeta*[mj(2:plength),mj(plength)];
    futureval0 = pbeta*[mj(1),mj(1:(plength-1))];
    
    futurevalderiv = (futureval1-futureval0)./(gap);
    futurevalderiv = repmat(futurevalderiv,plength,1);

    % horizontal: choices
    % vertical: state k 
    
    I = (vgridk-(1-pdelta)*k');
    vfoc = futurevalderiv - (p+p*pcvx*(I./k'));  

    % interpolation

    low = sum(vfoc>0,2);
    low(low<1) = 1;
    low(low>=(length(vfoc)-1)) = length(vfoc)-1;
    high = low + 1;
        
    indxlow     = sub2ind(size(vfoc),(1:plength)',low);
    indxhigh    = sub2ind(size(vfoc),(1:plength)',high);
    weightLow   = (vfoc(indxhigh) - 0)...
                ./(vfoc(indxhigh)-vfoc(indxlow));
    weightLow(weightLow<0) = 0;
    weightLow(weightLow>1) = 1;
    kprimeopt   = weightLow'.*vgridk(low)+(1-weightLow').*vgridk(high);    
    
    lb = min(vgridk);
    ub = max(vgridk);
    
    kprimeopt(kprimeopt<lb) = lb;
    kprimeopt(kprimeopt>ub) = ub;

end

  