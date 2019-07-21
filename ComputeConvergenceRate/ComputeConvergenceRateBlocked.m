function ConvergenceRate = ComputeConvergenceRateBlocked(MNmax,smax,ds)

% intialise
if nargin < 2, s = 2; end
if nargin < 3, ds = 0.05; end

s = smax:-ds:0; Ns = numel(s);
p = normcdf(s);

% compute roots
for m = 1:MNmax
    for n = 1:MNmax
        if m <= n % solve only if m<n
            P_down = [p(1:Ns-1).^n 0]/n;
            P_up = [0 sum(cell2mat(arrayfun(@(i) nchoosek(n,i)*((1-p(2:Ns)).^i).*p(2:Ns).^(n-i),m:n,'UniformOutput',false)'),1)]/n;
            P_stay = 1 - (P_down + P_up);
            TPM = gallery('tridiag',P_up(2:Ns),P_stay,P_down(1:Ns-1)); % TPM = transition probability matrix
            lambda = eig(full(TPM)); lambda = sort(abs(lambda),'descend'); % sorted eigenvals
            ConvergenceRate(m,n) = lambda(2);
        end
    end
end