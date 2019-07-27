classdef nonconsecutive < handle
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
        function this = nonconsecutive(MN_max,gamma,s_max,ds,responsemodel)
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
                    this.BalancePoint.Percentile(m,n) = (n*gamma)/(m + (n*gamma));
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
                    for k = 2:Ns, lambda(k) = (p(k-1)/(1-p(k)))*(m/n); end
                    lambda(1) = 1/(1 + sum(cumprod(lambda(2:end))));
                    if lambda(1)==0, lambda(1)=1e-300; end
                    this.StationaryDistribution.x_stim(m,n,:) = s; this.StationaryDistribution.x_prob(m,n,:) = p;
                    this.StationaryDistribution.y(m,n,:) = cumprod(lambda)/sum(cumprod(lambda));
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
                for n = 1:M_max
                    P_down = [p(1:Ns-1)/n 0];
                    P_up = [0 (1-p(2:Ns))/m];
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