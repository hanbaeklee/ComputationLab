function [kprimeopt] = fnoptinvirreversible(pbeta,vgridk,p,tempmat,mlambda)

    plength = length(vgridk);
    vfoc = pbeta*tempmat + mlambda - p;  

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

  