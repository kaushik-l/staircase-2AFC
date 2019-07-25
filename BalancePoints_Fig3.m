% load staircase solutions
load('staircases.mat');

%% define colors
clrs = [0 0.4 0.8 ; 0 0.8 0.4 ; 0.8 0 0.8 ; 0.8 0.4 0];
figure; hold on; set(gcf,'Position',[300 300 600 934]);

%% 1-up N-down
M = 1;
subplot(3,2,1); hold on; plot(1:NonconsecutiveStaircase.DesignParameter.N_max,...
    squeeze(NonconsecutiveStaircase.BalancePoint.Percentile(M,:)),'-','Color',clrs(1,:),'MarkerFaceColor',clrs(1,:)); 
subplot(3,2,1); hold on; plot(1:ConsecutiveStaircase.DesignParameter.N_max,...
    squeeze(ConsecutiveStaircase.BalancePoint.Percentile(M,:)),'-','Color',clrs(2,:),'MarkerFaceColor',clrs(2,:)); 
subplot(3,2,1); hold on; plot(1:BlockedStaircase.DesignParameter.N_max,...
    squeeze(BlockedStaircase.BalancePoint.Percentile(M,:)),'--','Color',clrs(3,:),'MarkerFaceColor',clrs(3,:)); 
subplot(3,2,1); hold on; plot(1:ContinuousStaircase.DesignParameter.N_max,...
    squeeze(ContinuousStaircase.BalancePoint.Percentile(M,:)),'-','Color',clrs(4,:),'MarkerFaceColor',clrs(4,:)); 
axis([1 NonconsecutiveStaircase.DesignParameter.N_max 0 1]);
set(gca,'YTick',0:0.25:1); hline(0.5,'--k');
set(gca,'TickDir','Out'); ticklength;
title('1-up N-down','fontweight','bold'); xlabel('\it N'); ylabel({'Cum. prob. of balance point, \Phi(s*)'}); 
legend('Nonconsecutive', 'Consecutive', 'Blocked', 'Continuous', 'Location','se'); hleg = legend('show'); title(hleg,'Update algorithm');

%% M-up 1-down
N = 1; 
subplot(3,2,2); hold on; plot(1:NonconsecutiveStaircase.DesignParameter.M_max,...
    squeeze(NonconsecutiveStaircase.BalancePoint.Percentile(:,N)),'-','Color',clrs(1,:),'MarkerFaceColor',clrs(1,:)); 
subplot(3,2,2); hold on; plot(1:ConsecutiveStaircase.DesignParameter.M_max,...
    squeeze(ConsecutiveStaircase.BalancePoint.Percentile(:,N)),'-','Color',clrs(2,:),'MarkerFaceColor',clrs(2,:)); 
subplot(3,2,2); hold on; plot(1:BlockedStaircase.DesignParameter.M_max,...
    squeeze(BlockedStaircase.BalancePoint.Percentile(:,N)),'--','Color',clrs(3,:),'MarkerFaceColor',clrs(3,:)); 
subplot(3,2,2); hold on; plot(1:ContinuousStaircase.DesignParameter.M_max,...
    squeeze(ContinuousStaircase.BalancePoint.Percentile(:,N)),'-','Color',clrs(4,:),'MarkerFaceColor',clrs(4,:)); 
axis([1 NonconsecutiveStaircase.DesignParameter.M_max 0 1]);
set(gca,'YTick',0:0.25:1); hline(0.5,'--k');
set(gca,'TickDir','Out'); ticklength;
title('M-up 1-down','fontweight','bold'); xlabel('\it M'); 


%% M-up N-down
clr_range = [0.05 0.95];
subplot(3,2,3); imagesc(1:NonconsecutiveStaircase.DesignParameter.M_max,1:NonconsecutiveStaircase.DesignParameter.N_max,...
    NonconsecutiveStaircase.BalancePoint.Percentile,clr_range); 
set(gca,'YDir','normal'); xlabel('\it N'); ylabel({'\it M'}); title('Nonconsecutive','fontweight','bold');
set(gca,'TickDir','Out'); ticklength;
subplot(3,2,4); imagesc(1:ConsecutiveStaircase.DesignParameter.M_max,1:ConsecutiveStaircase.DesignParameter.N_max,...
    ConsecutiveStaircase.BalancePoint.Percentile,clr_range); 
set(gca,'YDir','normal'); xlabel('\it N'); ylabel({'\it M'}); title('Consecutive','fontweight','bold');
set(gca,'TickDir','Out'); ticklength;
subplot(3,2,5); imagesc(1:BlockedStaircase.DesignParameter.M_max,1:BlockedStaircase.DesignParameter.N_max,...
    BlockedStaircase.BalancePoint.Percentile,clr_range); 
set(gca,'YDir','normal'); xlabel('\it N'); ylabel({'\it M'}); title('Blocked','fontweight','bold');
set(gca,'TickDir','Out'); ticklength;
cb = colorbar('northoutside','AxisLocation','in','Ticks',[0.25 0.5 0.75]); cb.Position = cb.Position + [0.27 0.08 -0.1 0]; 
cb.Label.String = '\Phi(s*)'; cb.Label.FontSize = 10; cb.Label.Position = cb.Label.Position + [0 4 0];
subplot(3,2,6); imagesc(1:ContinuousStaircase.DesignParameter.M_max,1:ContinuousStaircase.DesignParameter.N_max,...
    ContinuousStaircase.BalancePoint.Percentile,clr_range); 
set(gca,'YDir','normal'); xlabel('\it N'); ylabel({'\it M'}); title('Continuous','fontweight','bold');
set(gca,'TickDir','Out'); ticklength;
cmap = goodcolormap('bwr',10000); 
colormap(hot);