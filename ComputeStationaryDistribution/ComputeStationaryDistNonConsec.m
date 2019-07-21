function [StationaryDist,StationaryDistEntropy] = ComputeStationaryDistNonConsec(MNmax,smax,ds)

% intialise
if nargin < 2, s = 2; end
if nargin < 3, ds = 0.05; end

s = smax:-ds:-smax; Ns = numel(s);
p = normcdf(s);
StationaryDist = zeros(MNmax,MNmax,Ns);
StationaryDistEntropy = zeros(MNmax,MNmax);

% compute roots
for m = 1:MNmax
    lambda = nan(1,Ns);
    for n = 1:MNmax
        if m <= n % solve only if m<n
            for k = 2:Ns, lambda(k) = (p(k-1)/(1-p(k)))*(m/n); end
            lambda(1) = 1/(1 + sum(cumprod(lambda(2:end))));
            if lambda(1)==0, lambda(1)=1e-300; end
            StationaryDist(m,n,:) = cumprod(lambda)/sum(cumprod(lambda));
            % entropy of stationary distribution
            StatDist = squeeze(StationaryDist(m,n,:)); LogStatDist = log2(StatDist);
            StationaryDistEntropy(m,n) = -sum(StatDist(LogStatDist>-300).*LogStatDist(LogStatDist>-300));
        end
    end
end