function ConvergenceRate = ComputeConvergenceRateConsec(MNmax,smax,ds)

% intialise
if nargin < 2, s = 2; end
if nargin < 3, ds = 0.05; end

s = smax:-ds:0; Ns = numel(s);
p = normcdf(s);

% compute roots
for m = 1:MNmax
    for n = 1:MNmax
        if m <= n % solve only if m<n
            P_down = [((1-p(1:Ns-1)).*p(1:Ns-1).^n)./(1-p(1:Ns-1).^n) 0];
            P_up = [0 (p(2:Ns).*(1-p(2:Ns)).^m)./(1-(1-p(2:Ns)).^m)];
            P_stay = 1 - (P_down + P_up);
            TPM = gallery('tridiag',P_up(2:Ns),P_stay,P_down(1:Ns-1)); % TPM = transition probability matrix
            lambda = eig(full(TPM)); lambda = sort(abs(lambda),'descend'); % sorted eigenvals
            ConvergenceRate(m,n) = lambda(2);
        end
    end
end