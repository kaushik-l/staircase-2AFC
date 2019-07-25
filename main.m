%%
% N_objects = 1;
N_objects = 8;
MN_max = 10*ones(1,N_objects);
gamma = 1*ones(1,N_objects);
s_max = 2*ones(1,N_objects);
% ds = [0.1];
ds = [0.01 0.05 0.1 0.2:0.2:1];

%%
for k=1:numel(ds)
    NonconsecutiveStaircase(k) = nonconsecutive(MN_max(k),gamma(k),s_max(k),ds(k));
    NonconsecutiveStaircase(k).ComputeBalacePoint;
    NonconsecutiveStaircase(k).ComputeStationaryDistribution;
    NonconsecutiveStaircase(k).ComputeConvergenceRate;
    
    ConsecutiveStaircase(k) = consecutive(MN_max(k),gamma(k),s_max(k),ds(k));
    ConsecutiveStaircase(k).ComputeBalacePoint;
    ConsecutiveStaircase(k).ComputeStationaryDistribution;
    ConsecutiveStaircase(k).ComputeConvergenceRate;
    
    BlockedStaircase(k) = blocked(MN_max(k),gamma(k),s_max(k),ds(k));
    BlockedStaircase(k).ComputeBalacePoint;
    BlockedStaircase(k).ComputeStationaryDistribution;
    BlockedStaircase(k).ComputeConvergenceRate;
    
    ContinuousStaircase(k) = continuous(MN_max(k),gamma(k),s_max(k),ds(k));
    ContinuousStaircase(k).ComputeBalacePoint;
    ContinuousStaircase(k).ComputeStationaryDistribution;
    ContinuousStaircase(k).ComputeConvergenceRate;
end