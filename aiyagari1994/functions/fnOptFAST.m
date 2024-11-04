function [aprimeopt] = fnOptFAST(pBbeta,budget,vGrida,tempValue,pNumGrida,minVal)

% optimization option
options = optimset('Display','off','TolFun',1e-12,'TolX',1e-12,'MaxIter',20000,'MaxFunEvals',20000); % MaxFunEvals: 22400 (Default);
[aprimeopt] = fminbnd(@residCal,minVal,budget,options);

%=========================    
%nested function
%=========================    
function [optVal] = residCal(Guess)
    
aprime = Guess;
aLow = sum(vGrida<aprime);
aLow(aLow<1) = 1;
aLow(aLow>=pNumGrida) = pNumGrida-1;
aHigh = aLow+1;

weightLow = (vGrida(aHigh) - aprime)/(vGrida(aHigh)-vGrida(aLow));
weightLow(weightLow<0) = 0;
weightLow(weightLow>1) = 1;

value = weightLow*tempValue(aLow)+(1-weightLow)*tempValue(aHigh);
c = budget - aprime;

%update
value = log(c)+pBbeta*value;
optVal = -value;

end    
%=========================    
%=========================    

end