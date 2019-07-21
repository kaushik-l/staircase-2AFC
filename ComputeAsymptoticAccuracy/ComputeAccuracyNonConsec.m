function accuracy = ComputeAccuracyNonConsec(MNmax,gamma)

% intialise
if nargin < 2, gamma = 1; end
accuracy = zeros(MNmax,MNmax);

% compute
for m = 1:MNmax
    for n = 1:MNmax
        accuracy(m,n) = (n*gamma)/(m + (n*gamma));
    end
end