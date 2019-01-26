function accuracy = Compute_AccuracyContinuous(MNmax,plt)

% intialise
if nargin < 2, plt = false; end
accuracy = zeros(MNmax,MNmax);

% compute roots
for m = 1:MNmax
    for n = 1:MNmax
        if m < n % solve only if m<n
            rootvec = zeros(1,n+1); % Cofficients C of C(1)*p^n + C(2)*p^(n-1) + ... + C(n)*p + C(n+1)*1 = 0
            rootvec(end-n) = 1; % coefficient of p^n
            for i = 0:m
                rootvec(end-i) = rootvec(end-i) + (-1)^(i+1)*factorial(m)/(factorial(i)*factorial(m-i));
            end
            allroots = roots(rootvec);
            % pick only real roots between 0.5 and 1
            accuracy(m,n) = real(allroots(abs(imag(allroots)) < 1e-3 & allroots > 0.5 & allroots < 1));
        else, accuracy(m,n) = 0.5;  % accuracy is always 0.5 if m>=n
        end
    end
end

% plot
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
    % draw a triangle around the region m > n
    hold on;
    plot(1:MNmax,1:MNmax,'k')
    plot(1:MNmax,MNmax*ones(1,MNmax),'k')
    plot(ones(1,MNmax),1:MNmax,'k')
    text(1.5,0.7*MNmax,'Chance (50%)','fontsize',16)    
    % labels
    title('Accuracy Threshold (Continuous)','fontsize',20)
    xlabel('N','fontsize',14), ylabel('M','fontsize',14)
    set(gca,'fontsize',14)
    axis([0 MNmax+1 0 MNmax+1])
    axis square
end