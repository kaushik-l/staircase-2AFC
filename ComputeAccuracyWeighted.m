function accuracy = ComputeAccuracyWeighted(gammamax)

% compute roots
gamma = (1:gammamax);
accuracy = gamma./(1 + gamma);