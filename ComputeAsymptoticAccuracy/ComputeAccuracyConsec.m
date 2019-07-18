function accuracy = ComputeAccuracyConsec(MNmax,gamma)

% intialise
if nargin < 2, gamma = 1; end
accuracy = zeros(MNmax,MNmax);

% compute roots
for m = 1:MNmax
    for n = 1:MNmax
        if m < n
            rootvec = zeros(1,n+m); % Cofficients C of C(1)*p^(m+n-1) + C(2)*p^(m+n-2) + ... + C(m+n-1)*p + C(m+n)*1 = 0
            rootvec(end-(n-1)) = 1; % coefficient of p^(n-1)
            for i = 0:m
                % coefficients of p^(N-1) * (1-p)^M term
                rootvec(end-(i+n-1)) = rootvec(end-(i+n-1)) + (-1)^(i+1)*factorial(m)/(factorial(i)*factorial(m-i));
            end
            for i = 0:m-1
                % coefficients of (1-p)^(M-1) term
                rootvec(end-i) = rootvec(end-i) + gamma*(-1)^(i+1)*factorial(m-1)/(factorial(i)*factorial(m-1-i));
                % coefficients of p^N * (1-p)^(M-1) term
                rootvec(end-(i+n)) = rootvec(end-(i+n)) + gamma*(-1)^(i)*factorial(m-1)/(factorial(i)*factorial(m-1-i));
            end
            allroots = roots(rootvec);
            % pick only real roots between 0.5 and 1
            realroot = real(allroots(abs(imag(allroots)) < 1e-3 & real(allroots) >= 0.499 & real(allroots) < 1));
            if ~isempty(realroot), accuracy(m,n) = realroot; else, accuracy(m,n) = 0.5; end
        else, accuracy(m,n) = 0.5;   % accuracy is always 0.5 if m>=n
        end
    end
end