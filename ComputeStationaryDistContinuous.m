function StationaryDist = ComputeStationaryDistContinuous(MNmax,ds,gamma)

% intialise
if nargin < 2, ds = 0.1; end
if nargin < 3, gamma = 1; end

s = ds:ds:100*ds; Ns = numel(s);
Ps = erf(s./sqrt(2))/2 + 0.5;
StationaryDist = zeros(MNmax,MNmax,Ns);

% compute roots
for m = 1:MNmax
    for n = 1:MNmax
        if m < n % solve only if m<n
            rootvec = zeros(1,n+1); % Cofficients C of C(1)*p^n + C(2)*p^(n-1) + ... + C(n)*p + C(n+1)*1 = 0
            rootvec(end-n) = 1; % coefficient of p^n
            for i = 0:m
                rootvec(end-i) = rootvec(end-i) + gamma*(-1)^(i+1)*factorial(m)/(factorial(i)*factorial(m-i));
            end
            allroots = roots(rootvec);
            % pick only real roots between 0.5 and 1
            realroot = real(allroots(abs(imag(allroots)) < 1e-3 & real(allroots) >= 0.499 & real(allroots) < 1));
            if ~isempty(realroot), accuracy(m,n) = realroot; else, accuracy(m,n) = 0.5; end
        else, accuracy(m,n) = 0.5;  % accuracy is always 0.5 if m>=n
        end
    end
end