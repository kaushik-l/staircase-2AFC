function accuracy = ComputeAccuracyBlocked(MNmax,gamma,plt)

% intialise
if nargin < 2, gamma = 1; plt = false; end
if nargin < 3, plt = false; end
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

%plot
if plt
    figure; hold on;
    % draw circles proportional to accuracy
    for m = 1:MNmax
        for n = 1:MNmax
            if m < n
                plot(n,m,'ro','markersize',20)
                plot(n,m,'bo','markersize',20*(2*accuracy(m,n)-1),'markerfacecolor','b')
            end
        end
    end
    %
    hold on
    plot(1:MNmax,1:MNmax,'k')
    plot(1:MNmax,MNmax*ones(1,MNmax),'k')
    plot(ones(1,MNmax),1:MNmax,'k')
    text(1.5,0.7*MNmax,'Chance (50%)','fontsize',16)
    %
    title('Accuracy Threshold (Blocked)','fontsize',20)
    xlabel('N','fontsize',14), ylabel('M','fontsize',14)
    set(gca,'fontsize',14)
    axis([0 MNmax+1 0 MNmax+1])
    axis square
end