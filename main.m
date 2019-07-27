%%
clear;
MN_max = 10;
gamma = 1;
s_max = 2;
ds = 0.1;
responsemodel = 'gaussian';

NonconsecutiveStaircase = nonconsecutive(MN_max,gamma,s_max,ds,responsemodel);
NonconsecutiveStaircase.ComputeBalacePoint;
NonconsecutiveStaircase.ComputeStationaryDistribution;
NonconsecutiveStaircase.ComputeConvergenceRate;

ConsecutiveStaircase = consecutive(MN_max,gamma,s_max,ds,responsemodel);
ConsecutiveStaircase.ComputeBalacePoint;
ConsecutiveStaircase.ComputeStationaryDistribution;
ConsecutiveStaircase.ComputeConvergenceRate;

BlockedStaircase = blocked(MN_max,gamma,s_max,ds,responsemodel);
BlockedStaircase.ComputeBalacePoint;
BlockedStaircase.ComputeStationaryDistribution;
BlockedStaircase.ComputeConvergenceRate;

ContinuousStaircase = continuous(MN_max,gamma,s_max,ds,responsemodel);
ContinuousStaircase.ComputeBalacePoint;
ContinuousStaircase.ComputeStationaryDistribution;
ContinuousStaircase.ComputeConvergenceRate;

save('staircases.mat', 'NonconsecutiveStaircase','ConsecutiveStaircase','BlockedStaircase','ContinuousStaircase');

%%
clear;
MN_max = 10;
gamma = 1;
s_max = 2;
ds = [0.01 0.05 0.1 0.2:0.2:1];
responsemodel = 'gaussian';

for k=1:numel(ds)
    NonconsecutiveStaircase(k) = nonconsecutive(MN_max,gamma,s_max,ds(k),responsemodel);
    NonconsecutiveStaircase(k).ComputeBalacePoint;
    NonconsecutiveStaircase(k).ComputeStationaryDistribution;
    NonconsecutiveStaircase(k).ComputeConvergenceRate;
    
    ConsecutiveStaircase(k) = consecutive(MN_max,gamma,s_max,ds(k),responsemodel);
    ConsecutiveStaircase(k).ComputeBalacePoint;
    ConsecutiveStaircase(k).ComputeStationaryDistribution;
    ConsecutiveStaircase(k).ComputeConvergenceRate;
    
    BlockedStaircase(k) = blocked(MN_max,gamma,s_max,ds(k),responsemodel);
    BlockedStaircase(k).ComputeBalacePoint;
    BlockedStaircase(k).ComputeStationaryDistribution;
    BlockedStaircase(k).ComputeConvergenceRate;
    
    ContinuousStaircase(k) = continuous(MN_max,gamma,s_max,ds(k),responsemodel);
    ContinuousStaircase(k).ComputeBalacePoint;
    ContinuousStaircase(k).ComputeStationaryDistribution;
    ContinuousStaircase(k).ComputeConvergenceRate;
end

save('staircases_vs_stepsize.mat', 'NonconsecutiveStaircase','ConsecutiveStaircase','BlockedStaircase','ContinuousStaircase');

%%
clear;
MN_max = 10;
gamma = 1;
s_max = [2 2 4 2];
ds = 0.1;
responsemodel = {'gaussian','logistic','weibull','linear'};

for k=1:numel(responsemodel)
    NonconsecutiveStaircase(k) = nonconsecutive(MN_max,gamma,s_max(k),ds,responsemodel{k});
    NonconsecutiveStaircase(k).ComputeBalacePoint;
    NonconsecutiveStaircase(k).ComputeStationaryDistribution;
    NonconsecutiveStaircase(k).ComputeConvergenceRate;
    
    ConsecutiveStaircase(k) = consecutive(MN_max,gamma,s_max(k),ds,responsemodel{k});
    ConsecutiveStaircase(k).ComputeBalacePoint;
    ConsecutiveStaircase(k).ComputeStationaryDistribution;
    ConsecutiveStaircase(k).ComputeConvergenceRate;
    
    BlockedStaircase(k) = blocked(MN_max,gamma,s_max(k),ds,responsemodel{k});
    BlockedStaircase(k).ComputeBalacePoint;
    BlockedStaircase(k).ComputeStationaryDistribution;
    BlockedStaircase(k).ComputeConvergenceRate;
    
    ContinuousStaircase(k) = continuous(MN_max,gamma,s_max(k),ds,responsemodel{k});
    ContinuousStaircase(k).ComputeBalacePoint;
    ContinuousStaircase(k).ComputeStationaryDistribution;
    ContinuousStaircase(k).ComputeConvergenceRate;
end

save('staircases_vs_responsemodel.mat', 'NonconsecutiveStaircase','ConsecutiveStaircase','BlockedStaircase','ContinuousStaircase');