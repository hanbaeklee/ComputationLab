function output = fnSimulator(iniPoint,mTrans,simTime)

cum   = cumsum(mTrans')';
output    = ones(simTime,1);
output(1) = iniPoint;

for it = 2:simTime
    x = find(rand <= cum(output(it-1),:));
    output(it) = x(1);
end
