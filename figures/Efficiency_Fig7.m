clear;
load staircases_vs_responsemodel.mat
nummodels = numel(NonconsecutiveStaircase);

%% gaussian/logistic/weibull/linear
% define colors
clrs = [0 0.4 0.8 ; 0 0.8 0.4 ; 0.8 0 0.8 ; 0.8 0.4 0];

for k=1:nummodels
    figure; hold on; set(gcf,'Position',[50 300 1000 1000]);

    % 1-up 1-down, 10-up 10-down
    M = 1; N = 1;
    subplot(3,3,1); hold on;
    plot(squeeze(ConsecutiveStaircase(k).StationaryDistribution.x_stim(M,N,:)),...
        squeeze(ConsecutiveStaircase(k).StationaryDistribution.y(M,N,:)),'Color',clrs(2,:));
    set(gca,'XTick',[-2 0 2],'XTickLabel',{'-2\sigma','0','+2\sigma'},'YTick',[0 0.25 0.5]);
    set(gca,'TickDir','Out'); ticklength;
    xlabel('Stimulus, \it s'); ylabel('Stationary probability, \pi(s)'); title('Stationary distribution','Fontweight','Bold');
    M = 10; N = 10;
    plot(squeeze(ConsecutiveStaircase(k).StationaryDistribution.x_stim(M,N,:)),...
        squeeze(ConsecutiveStaircase(k).StationaryDistribution.y(M,N,:)),'--','Color',clrs(2,:));
    legend('1-up 1-down','10-up 10-down','Location','best'); hleg = legend('show'); title(hleg,'Design parameters');

    % Information gain vs (M,N)
    subplot(3,3,2); hold on;
    plot(1:NonconsecutiveStaircase(k).DesignParameter.N_max,...
        diag(NonconsecutiveStaircase(k).StationaryDistribution.EntropyUniform) - diag(NonconsecutiveStaircase(k).StationaryDistribution.Entropy),'-',...
        'Color',clrs(1,:));
    plot(1:ConsecutiveStaircase(k).DesignParameter.N_max,...
        diag(ConsecutiveStaircase(k).StationaryDistribution.EntropyUniform) - diag(ConsecutiveStaircase(k).StationaryDistribution.Entropy),'-',...
        'Color',clrs(2,:));
    plot(1:BlockedStaircase(k).DesignParameter.N_max,...
        diag(BlockedStaircase(k).StationaryDistribution.EntropyUniform) - diag(BlockedStaircase(k).StationaryDistribution.Entropy),'--',...
        'Color',clrs(3,:));
    plot(1:ContinuousStaircase(k).DesignParameter.N_max,...
        diag(ContinuousStaircase(k).StationaryDistribution.EntropyUniform) - diag(ContinuousStaircase(k).StationaryDistribution.Entropy),'-',...
        'Color',clrs(4,:));
    axis([1 10 0 5.35]);
    legend('Nonconsecutive', 'Consecutive', 'Blocked', 'Continuous', 'Location','nw'); hleg = legend('show'); title(hleg,'Update rule');
    xlabel('\it M,N'); ylabel('Information gain [bits]');
    set(gca,'TickDir','Out'); ticklength;

    % Information gain vs M,N
    clr_range = [1.8 3.6];
    subplot(3,3,4); hold on;
    imagesc(1:NonconsecutiveStaircase(k).DesignParameter.M_max,1:NonconsecutiveStaircase(k).DesignParameter.N_max,...
        NonconsecutiveStaircase(k).StationaryDistribution.EntropyUniform - NonconsecutiveStaircase(k).StationaryDistribution.Entropy,clr_range);
    axis([0.5 10.5 0.5 10.5]); set(gca,'YDir','normal'); xlabel('\it M'); ylabel({'\it N'}); title('Nonconsecutive','fontweight','bold');
    set(gca,'TickDir','Out'); ticklength;
    subplot(3,3,5); hold on;
    imagesc(1:ConsecutiveStaircase(k).DesignParameter.M_max,1:ConsecutiveStaircase(k).DesignParameter.N_max,...
        ConsecutiveStaircase(k).StationaryDistribution.EntropyUniform - ConsecutiveStaircase(k).StationaryDistribution.Entropy,clr_range);
    axis([0.5 10.5 0.5 10.5]); set(gca,'YDir','normal'); xlabel('\it M'); ylabel({'\it N'}); title('Consecutive','fontweight','bold');
    set(gca,'TickDir','Out'); ticklength;
    subplot(3,3,7); hold on;
    imagesc(1:BlockedStaircase(k).DesignParameter.M_max,1:BlockedStaircase(k).DesignParameter.N_max,...
        BlockedStaircase(k).StationaryDistribution.EntropyUniform - BlockedStaircase(k).StationaryDistribution.Entropy,clr_range);
    axis([0.5 10.5 0.5 10.5]); set(gca,'YDir','normal'); xlabel('\it M'); ylabel({'\it N'}); title('Blocked','fontweight','bold');
    cb = colorbar('northoutside','AxisLocation','in','Ticks',[1.8 2.7 3.6]); cb.Position = cb.Position + [0.19 0.075 -0.1 0];
    cb.Label.String = 'Info. gain [bits]'; cb.Label.FontSize = 10; cb.Label.Position = cb.Label.Position + [0 4 0];
    set(gca,'TickDir','Out'); ticklength;
    subplot(3,3,8); hold on;
    imagesc(1:ContinuousStaircase(k).DesignParameter.M_max,1:ContinuousStaircase(k).DesignParameter.N_max,...
        ContinuousStaircase(k).StationaryDistribution.EntropyUniform - ContinuousStaircase(k).StationaryDistribution.Entropy,clr_range);
    axis([0.5 10.5 0.5 10.5]); set(gca,'YDir','normal'); xlabel('\it M'); ylabel({'\it N'}); title('Continuous','fontweight','bold');
    set(gca,'TickDir','Out'); ticklength;
    cmap = goodcolormap('bwr',10000);
    colormap(hot);
end

%% summary
clrs = [0 0.4 0.8 ; 0 0.8 0.4 ; 0.8 0 0.8 ; 0.8 0.4 0];
figure; hold on; set(gcf,'Position',[50 300 800 500]);

for i=1:nummodels
    NumTrials = NonconsecutiveStaircase(i).ConvergenceRate.numtrials;
    InfoGain = NonconsecutiveStaircase(i).StationaryDistribution.EntropyUniform - NonconsecutiveStaircase(i).StationaryDistribution.Entropy;
    Efficiency(i,1,:) = InfoGain(:)./NumTrials(:);
    
    NumTrials = ConsecutiveStaircase(i).ConvergenceRate.numtrials;
    InfoGain = ConsecutiveStaircase(i).StationaryDistribution.EntropyUniform - ConsecutiveStaircase(i).StationaryDistribution.Entropy;
    Efficiency(i,2,:) = InfoGain(:)./NumTrials(:);
    
    NumTrials = BlockedStaircase(i).ConvergenceRate.numtrials;
    InfoGain = BlockedStaircase(i).StationaryDistribution.EntropyUniform - BlockedStaircase(i).StationaryDistribution.Entropy;
    Efficiency(i,3,:) = InfoGain(:)./NumTrials(:);
    
    NumTrials = ContinuousStaircase(i).ConvergenceRate.numtrials;
    InfoGain = ContinuousStaircase(i).StationaryDistribution.EntropyUniform - ContinuousStaircase(i).StationaryDistribution.Entropy;
    Efficiency(i,4,:) = InfoGain(:)./NumTrials(:);
    
    MeanEfficiency = squeeze(mean(Efficiency,3));
    SEMEfficiency = squeeze(std(Efficiency,[],3)/sqrt(100));    
end

axis off;
axes('Position',[0.13 0.3 0.78 0.6]); hold on;
h = bar(1:4,MeanEfficiency,0.7); 
for k=1:4, h(k).FaceColor = clrs(k,:); end
[ngroups,nbars] = size(MeanEfficiency);
BarWidth = h(1).BarWidth;
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, MeanEfficiency(:,i), SEMEfficiency(:,i), '.k');
    er.CapSize = 0; 
end
ylim([0 max(MeanEfficiency(:)) + max(SEMEfficiency(:))]); set(gca,'XTick',[],'YTick',0:0.01:max(MeanEfficiency(:)) + max(SEMEfficiency(:)));
ylabel('Average efficiency [bits / trial]');
legend('Nonconsecutive', 'Consecutive', 'Blocked', 'Continuous', 'Location','nw','Color','None'); hleg = legend('show'); title(hleg,'Update algorithm');
set(gca,'Color','none');

axes('Position',[0.165 0.1 0.12 0.12]); hold on; 
plot(squeeze(NonconsecutiveStaircase(1).StationaryDistribution.x_stim(1,1,:)),...
    squeeze(NonconsecutiveStaircase(1).StationaryDistribution.x_prob(1,1,:))); 
set(gca,'Color','none','YColor','none','XTick',[]); vline(0,'-k');
title('\rm Cumulative Gaussian','Fontsize',10);

axes('Position',[0.365 0.1 0.12 0.12]); hold on;
plot(squeeze(NonconsecutiveStaircase(2).StationaryDistribution.x_stim(1,1,:)),...
    squeeze(NonconsecutiveStaircase(2).StationaryDistribution.x_prob(1,1,:))); 
set(gca,'Color','none','YColor','none','XTick',[]); vline(0,'-k');
title('\rm Logistic','Fontsize',10);

axes('Position',[0.56 0.1 0.12 0.12]); hold on;
plot(squeeze(NonconsecutiveStaircase(3).StationaryDistribution.x_stim(1,1,:)),...
    squeeze(NonconsecutiveStaircase(3).StationaryDistribution.x_prob(1,1,:))); 
set(gca,'Color','none','XTick',[],'YTick',[]);
title('\rm Weibull','Fontsize',10);

axes('Position',[0.76 0.1 0.12 0.12]); hold on;
plot(squeeze(NonconsecutiveStaircase(4).StationaryDistribution.x_stim(1,1,:)),...
    squeeze(NonconsecutiveStaircase(4).StationaryDistribution.x_prob(1,1,:))); 
set(gca,'Color','none','YColor','none','XTick',[]); vline(0,'-k');
title('\rm Linear','Fontsize',10);