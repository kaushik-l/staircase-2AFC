function accuracy = compareStairThreshSolutionBlocked(MNmax,plt)

% intialise
if nargin < 2, plt = false; end
accuracy = zeros(MNmax,MNmax);

% compute roots
for M = 1:MNmax
    for N = 1:MNmax
        if M < N
            rootvec = zeros(1,N+M); % Cofficients C of C(1)*p^(m+n-1) + C(2)*p^(m+n-2) + ... + C(m+n-1)*p + C(m+n)*1 = 0
            rootvec(end-(N-1)) = 1; % coefficient of p^(n-1)
            for i = 0:M
                % coefficients of p^(N-1) * (1-p)^M term
                rootvec(end-(i+N-1)) = rootvec(end-(i+N-1)) + (-1)^(i+1)*factorial(M)/(factorial(i)*factorial(M-i));
            end
            for i = 0:M-1
                % coefficients of (1-p)^(M-1) term
                rootvec(end-i) = rootvec(end-i) + (-1)^(i+1)*factorial(M-1)/(factorial(i)*factorial(M-1-i));
                % coefficients of p^N * (1-p)^(M-1) term
                rootvec(end-(i+N)) = rootvec(end-(i+N)) + (-1)^(i)*factorial(M-1)/(factorial(i)*factorial(M-1-i));
            end
            allroots = roots(rootvec);
            % pick only real roots between 0.5 and 1
            accuracy(M,N) = real(allroots(abs(imag(allroots)) < 1e-3 & allroots > 0.5 & allroots < 1));
        else, accuracy(M,N) = 0.5;   % accuracy is always 0.5 if m>=n
        end
    end
end

if plt
    figure; hold on;
    %
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