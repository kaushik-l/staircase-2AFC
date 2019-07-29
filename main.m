%% basic sript for a pre-defined step size and response model
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

%% use different step sizes
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

%% use different response models
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

%% Simulate staircase
clear;
M = 1;
N = 2;
gamma = 1;
s_max = 2;
ds = 0.1;
responsemodel = 'gaussian';
numtrials = 10000;
numsimulations = 100;

% simulated
for k=1:numsimulations
    fprintf(['Stimulation # ' num2str(k) '\n']);
    
    NonconsecutiveStaircase(k) = nonconsecutive(max(M,N),gamma,s_max,ds,'gaussian');
    NonconsecutiveStaircase(k).SimulateProcess(M,N,numtrials);
    Simulated.NonconsecutiveStaircase.BalancePoint(k) = NonconsecutiveStaircase(k).Simulation.BalancePoint;
    
    ConsecutiveStaircase(k) = consecutive(max(M,N),gamma,s_max,ds,'gaussian');
    ConsecutiveStaircase(k).SimulateProcess(M,N,numtrials);
    Simulated.ConsecutiveStaircase.BalancePoint(k) = ConsecutiveStaircase(k).Simulation.BalancePoint;
    
    BlockedStaircase(k) = blocked(max(M,N),gamma,s_max,ds,'gaussian');
    BlockedStaircase(k).SimulateProcess(M,N,numtrials);
    Simulated.BlockedStaircase.BalancePoint(k) = BlockedStaircase(k).Simulation.BalancePoint;
    
    ContinuousStaircase(k) = continuous(max(M,N),gamma,s_max,ds,'gaussian');
    ContinuousStaircase(k).SimulateProcess(M,N,numtrials);
    Simulated.ContinuousStaircase.BalancePoint(k) = ContinuousStaircase(k).Simulation.BalancePoint;
end

% theoretical
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