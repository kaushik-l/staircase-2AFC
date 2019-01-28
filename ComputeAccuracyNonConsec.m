function accuracy = ComputeAccuracyNonConsec(MNmax,gamma,plt)

% intialise
if nargin < 2, gamma = 1; plt = false; end
if nargin < 3, plt = false; end
accuracy = zeros(MNmax,MNmax);

% compute
for m = 1:MNmax
    for n = 1:MNmax
        if m < n, accuracy(m,n) = (n*gamma)/(m + (n*gamma));
        else, accuracy(m,n) = 0.5; end
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
    %
    hold on
    plot(1:MNmax,1:MNmax,'k')
    plot(1:MNmax,MNmax*ones(1,MNmax),'k')
    plot(ones(1,MNmax),1:MNmax,'k')
    text(1.5,0.7*MNmax,'Chance (50%)','fontsize',16)
    %
    title('Accuracy Threshold (Non-consecutive)','fontsize',20)
    xlabel('N','fontsize',14), ylabel('M','fontsize',14)
    set(gca,'fontsize',14)
    axis([0 MNmax+1 0 MNmax+1])
    axis square    
end