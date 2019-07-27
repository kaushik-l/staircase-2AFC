classdef blocked < handle
    %%
    properties
        DesignParameter
        BalancePoint
        StationaryDistribution
        ConvergenceRate
    end
    %%
    methods
        
        %% class constructor
        function this = blocked(MN_max,gamma,s_max,ds,responsemodel)
            this.DesignParameter.M_max = MN_max;
            this.DesignParameter.N_max = MN_max;
            this.DesignParameter.gamma = gamma;
            this.DesignParameter.stim_max = s_max;
            this.DesignParameter.stepsize = ds;
            this.DesignParameter.responsemodel = responsemodel;
        end
        
        %% compute balance point
        function this = ComputeBalacePoint(this)
            
            M_max = this.DesignParameter.M_max; N_max = this.DesignParameter.N_max; gamma = this.DesignParameter.gamma;
            this.BalancePoint.Percentile = zeros(M_max,N_max);
            % compute roots
            for m = 1:M_max
                for n = 1:N_max
                    if m <= n % solve only if m<n
                        rootvec = zeros(1,n+1); % Cofficients C of C(1)*p^n + C(2)*p^(n-1) + ... + C(n)*p + C(n+1)*1 = 0
                        rootvec(end-n) = 1; % coefficient of p^n
                        for i = m:n
                            for j = 0:i
                                rootvec(end-(n-i+j)) = rootvec(end-(n-i+j)) + gamma*(-1)^(j+1)*factorial(n)/(factorial(n-i)*factorial(j)*factorial(i-j));
                            end
                        end
                    else
                        rootvec = zeros(1,m+1); % Cofficients C of C(1)*p^n + C(2)*p^(n-1) + ... + C(n)*p + C(n+1)*1 = 0
                        rootvec(end-m) = 1; % coefficient of p^m
                        for i = n:m
                            for j = 0:i
                                rootvec(end-(m-i+j)) = rootvec(end-(m-i+j)) + gamma*(-1)^(j+1)*factorial(m)/(factorial(m-i)*factorial(j)*factorial(i-j));
                            end
                        end
                    end
                    allroots = roots(rootvec);
                    % pick only real roots between 0 and 1
                    realroot = real(allroots(abs(imag(allroots)) < 1e-3 & real(allroots) >= 1e-3 & real(allroots) < 1));
                    if ~isempty(realroot), this.BalancePoint.Percentile(m,n) = realroot; else, this.BalancePoint.Percentile(m,n) = 0.5; end
                    if m>n, this.BalancePoint.Percentile(m,n) = 1-this.BalancePoint.Percentile(m,n); end %due to symmetry of the problem
                end
            end
            
        end
        
        %% compute stationary distribution
        function this = ComputeStationaryDistribution(this)
            
            M_max = this.DesignParameter.M_max; N_max = this.DesignParameter.N_max;
            smax = this.DesignParameter.stim_max; ds = this.DesignParameter.stepsize;
            switch lower(this.DesignParameter.responsemodel)
                case 'gaussian'
                    s = smax:-ds:-smax; % stimulus levels
                    Ns = numel(s); % number of stimulus levels
                    p = normcdf(s); % percentiles corresponding to stimulus levels
                case 'logistic'
                    s = smax:-ds:-smax; % stimulus levels
                    Ns = numel(s); % number of stimulus levels
                    p = (1./(1 + exp(-s*pi/sqrt(3)))); % percentiles corresponding to stimulus levels
                case 'weibull'
                    s = ds + (smax:-ds:0); % stimulus levels
                    Ns = numel(s); % number of stimulus levels
                    p = (1 - exp(-s)); % percentiles corresponding to stimulus levels
                case 'linear'
                    s = smax:-ds:-smax; % stimulus levels
                    Ns = numel(s); % number of stimulus levels
                    p = flip(cumsum(ones(1,Ns)/(Ns+1)));
            end
            this.StationaryDistribution.x_stim = zeros(M_max,N_max,Ns);
            this.StationaryDistribution.x_prob = zeros(M_max,N_max,Ns);
            this.StationaryDistribution.y = zeros(M_max,N_max,Ns);
            this.StationaryDistribution.Entropy = zeros(M_max,N_max);            
            % compute roots
            for m = 1:M_max
                lambda = nan(1,Ns);
                for n = 1:N_max
                    if m <= n % solve only if m<n
                        for k = 2:Ns, lambda(k) = p(k-1)^n/sum(arrayfun(@(i) nchoosek(n,i)*((1-p(k))^i).*p(k)^(n-i),m:n)); end
                    else
                        for k = 2:Ns, lambda(k) = p(k-1)^m/sum(arrayfun(@(i) nchoosek(m,i)*((1-p(k))^i).*p(k)^(m-i),n:m)); end
                    end
                    lambda(1) = 1/(1 + sum(cumprod(lambda(2:end))));
                    if lambda(1)==0, lambda(1)=1e-300; end
                    this.StationaryDistribution.y(m,n,:) = cumprod(lambda)/sum(cumprod(lambda));
                    if m<=n, this.StationaryDistribution.x_stim(m,n,:) = s; this.StationaryDistribution.x_prob(m,n,:) = p;
                    else, this.StationaryDistribution.x_stim(m,n,:) = -s; this.StationaryDistribution.x_prob(m,n,:) = 1-p; end
                    % entropy of stationary distribution
                    StatDist = squeeze(this.StationaryDistribution.y(m,n,:)); LogStatDist = log2(StatDist);
                    this.StationaryDistribution.Entropy(m,n) = -sum(StatDist(LogStatDist>-300).*LogStatDist(LogStatDist>-300));
                    this.StationaryDistribution.EntropyUniform(m,n) = -log2(1/Ns);
                end
            end
            
        end
        
        %% compute convergence rate
        function this = ComputeConvergenceRate(this)
            
            M_max = this.DesignParameter.M_max; N_max = this.DesignParameter.N_max;
            smax = this.DesignParameter.stim_max; ds = this.DesignParameter.stepsize;
            switch lower(this.DesignParameter.responsemodel)
                case 'gaussian'
                    s = smax:-ds:-smax; % stimulus levels
                    Ns = numel(s); % number of stimulus levels
                    p = normcdf(s); % percentiles corresponding to stimulus levels
                case 'logistic'
                    s = smax:-ds:-smax; % stimulus levels
                    Ns = numel(s); % number of stimulus levels
                    p = (1./(1 + exp(-s*pi/sqrt(3)))); % percentiles corresponding to stimulus levels
                case 'weibull'
                    s = ds + (smax:-ds:0); % stimulus levels
                    Ns = numel(s); % number of stimulus levels
                    p = (1 - exp(-s)); % percentiles corresponding to stimulus levels
                case 'linear'
                    s = smax:-ds:-smax; % stimulus levels
                    Ns = numel(s); % number of stimulus levels
                    p = flip(cumsum(ones(1,Ns)/(Ns+1)));
            end
            this.ConvergenceRate.x_stim = zeros(M_max,N_max,Ns);
            this.ConvergenceRate.x_prob = zeros(M_max,N_max,Ns);
            % compute roots
            for m = 1:M_max
                for n = 1:N_max
                    if m <= n % solve only if m<n
                        P_down = [p(1:Ns-1).^n 0]/n;
                        P_up = [0 sum(cell2mat(arrayfun(@(i) nchoosek(n,i)*((1-p(2:Ns)).^i).*p(2:Ns).^(n-i),m:n,'UniformOutput',false)'),1)]/n;
                    else
                        P_down = [p(1:Ns-1).^m 0]/m;
                        P_up = [0 sum(cell2mat(arrayfun(@(i) nchoosek(m,i)*((1-p(2:Ns)).^i).*p(2:Ns).^(m-i),n:m,'UniformOutput',false)'),1)]/m;
                    end
                    P_stay = 1 - (P_down + P_up);
                    TPM = gallery('tridiag',P_up(2:Ns),P_stay,P_down(1:Ns-1)); % TPM = transition probability matrix
                    lambda = abs(eig(full(TPM))); lambda(lambda>1-eps) = 1; lambda = unique(lambda);
                    lambda = sort(lambda,'descend'); % sorted eigenvals
                    if m==1 && n==1, this.ConvergenceRate.lambda(m,n) = lambda(3);
                    else, this.ConvergenceRate.lambda(m,n) = lambda(2); end
                    this.ConvergenceRate.numtrials(m,n) = ceil(log(1 - 0.99)/log(this.ConvergenceRate.lambda(m,n))); %99 percent convergence
                    this.ConvergenceRate.TPM(m,n,:,:) = full(TPM);
                    this.ConvergenceRate.x_stim(m,n,:) = s; this.ConvergenceRate.x_prob(m,n,:) = p;
                end
            end
            
        end            
    end
end