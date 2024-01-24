function [kprimeopt] = fnOptInv_matrix(pbeta,pdelta,pgamma_i,vgridk,p,tj0)

    [plength1, plength2] = size(tj0);
    vgridkorigup = [vgridk(2:plength1),vgridk(plength1)]';
    vgridkorigdn = [vgridk(1),vgridk(1:(plength1-1))]';
    gap = vgridkorigup-vgridkorigdn;
    gap = repmat(gap,1,plength2);
    
    futureval1 = pbeta*[tj0(2:plength1,:);tj0(plength1,:)];
    futureval0 = pbeta*[tj0(1,:);tj0(1:(plength1-1),:)];
    
    futurevalderiv = (futureval1-futureval0)./(gap);
    vfoc = futurevalderiv - p*pgamma_i;
    vfoc_trans = vfoc';

    % interpolation
    low = sum(vfoc>0,1);
    low(low<1) = 1;
    low(low>=(plength1-1)) = plength1-1;
    high = low+1;
    
    indxlow     = sub2ind(size(vfoc_trans),(1:plength2)',low');
    indxhigh    = sub2ind(size(vfoc_trans),(1:plength2)',high');
    weightLow   = (vfoc_trans(indxhigh) - 0)...
                ./(vfoc_trans(indxhigh)-vfoc_trans(indxlow));
    weightLow(weightLow<0) = 0;
    weightLow(weightLow>1) = 1;
    kprimeopt   = weightLow'.*vgridk(low)+(1-weightLow').*vgridk(low+1);    
    
    lb = min(vgridk);
    ub = max(vgridk);
    
    kprimeopt(kprimeopt<lb) = lb;
    kprimeopt(kprimeopt>ub) = ub;
    
end