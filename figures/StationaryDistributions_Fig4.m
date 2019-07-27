% load staircase solutions
clear;
load('staircases.mat');

%% define colors
clrs = [0 0.4 0.8 ; 0 0.8 0.4 ; 0.8 0 0.8 ; 0.8 0.4 0];
figure; hold on; set(gcf,'Position',[50 300 1000 1000]);

%% 1-up 1-down, 10-up 10-down
M = 1; N = 1;
subplot(3,3,1); hold on;
plot(squeeze(ConsecutiveStaircase.StationaryDistribution.x_stim(M,N,:)),...
    squeeze(ConsecutiveStaircase.StationaryDistribution.y(M,N,:)),'Color',clrs(2,:));
set(gca,'XTick',[-2 0 2],'XTickLabel',{'-2\sigma','0','+2\sigma'},'YTick',[0 0.25 0.5]);
set(gca,'TickDir','Out'); ticklength;
xlabel('Stimulus, \it s'); ylabel('Stationary probability, \pi(s)'); title('Stationary distribution','Fontweight','Bold');
M = 10; N = 10;
plot(squeeze(ConsecutiveStaircase.StationaryDistribution.x_stim(M,N,:)),...
    squeeze(ConsecutiveStaircase.StationaryDistribution.y(M,N,:)),'--','Color',clrs(2,:));
legend('1-up 1-down','10-up 10-down','Location','best','Color','none'); hleg = legend('show'); title(hleg,'Design parameters');
set(gca,'Color','none');

%% Information gain vs (M,N)
subplot(3,3,2); hold on;
plot(1:NonconsecutiveStaircase.DesignParameter.N_max,...
    diag(NonconsecutiveStaircase.StationaryDistribution.EntropyUniform) - diag(NonconsecutiveStaircase.StationaryDistribution.Entropy),'-',...
    'Color',clrs(1,:));
plot(1:ConsecutiveStaircase.DesignParameter.N_max,...
    diag(ConsecutiveStaircase.StationaryDistribution.EntropyUniform) - diag(ConsecutiveStaircase.StationaryDistribution.Entropy),'-',...
    'Color',clrs(2,:));
plot(1:BlockedStaircase.DesignParameter.N_max,...
    diag(BlockedStaircase.StationaryDistribution.EntropyUniform) - diag(BlockedStaircase.StationaryDistribution.Entropy),'--',...
    'Color',clrs(3,:));
plot(1:ContinuousStaircase.DesignParameter.N_max,...
    diag(ContinuousStaircase.StationaryDistribution.EntropyUniform) - diag(ContinuousStaircase.StationaryDistribution.Entropy),'-',...
    'Color',clrs(4,:));
axis([1 10 0 5.35]);
legend('Nonconsecutive', 'Consecutive', 'Blocked', 'Continuous', 'Location','nw','Color','none'); hleg = legend('show'); title(hleg,'Update rule');
xlabel('\it M,N'); ylabel('Information gain [bits]');
set(gca,'TickDir','Out'); ticklength;
set(gca,'Color','none');

%% Information gain vs M,N
clr_range = [1.8 3.6];
subplot(3,3,4); hold on;
imagesc(1:NonconsecutiveStaircase.DesignParameter.M_max,1:NonconsecutiveStaircase.DesignParameter.N_max,...
    NonconsecutiveStaircase.StationaryDistribution.EntropyUniform - NonconsecutiveStaircase.StationaryDistribution.Entropy,clr_range);
axis([0.5 10.5 0.5 10.5]); set(gca,'YDir','normal'); xlabel('\it M'); ylabel({'\it N'}); title('Nonconsecutive','fontweight','bold');
set(gca,'TickDir','Out'); ticklength; set(gca,'Color','none');
subplot(3,3,5); hold on;
imagesc(1:ConsecutiveStaircase.DesignParameter.M_max,1:ConsecutiveStaircase.DesignParameter.N_max,...
    ConsecutiveStaircase.StationaryDistribution.EntropyUniform - ConsecutiveStaircase.StationaryDistribution.Entropy,clr_range);
axis([0.5 10.5 0.5 10.5]); set(gca,'YDir','normal'); xlabel('\it M'); ylabel({'\it N'}); title('Consecutive','fontweight','bold');
set(gca,'TickDir','Out'); ticklength; set(gca,'Color','none');
subplot(3,3,7); hold on;
imagesc(1:BlockedStaircase.DesignParameter.M_max,1:BlockedStaircase.DesignParameter.N_max,...
    BlockedStaircase.StationaryDistribution.EntropyUniform - BlockedStaircase.StationaryDistribution.Entropy,clr_range);
axis([0.5 10.5 0.5 10.5]); set(gca,'YDir','normal'); xlabel('\it M'); ylabel({'\it N'}); title('Blocked','fontweight','bold');
cb = colorbar('northoutside','AxisLocation','in','Ticks',[1.8 2.7 3.6]); cb.Position = cb.Position + [0.19 0.075 -0.1 0]; 
cb.Label.String = 'Info. gain [bits]'; cb.Label.FontSize = 10; cb.Label.Position = cb.Label.Position + [0 4 0];
set(gca,'TickDir','Out'); ticklength; set(gca,'Color','none');
subplot(3,3,8); hold on;
imagesc(1:ContinuousStaircase.DesignParameter.M_max,1:ContinuousStaircase.DesignParameter.N_max,...
    ContinuousStaircase.StationaryDistribution.EntropyUniform - ContinuousStaircase.StationaryDistribution.Entropy,clr_range);
axis([0.5 10.5 0.5 10.5]); set(gca,'YDir','normal'); xlabel('\it M'); ylabel({'\it N'}); title('Continuous','fontweight','bold');
set(gca,'TickDir','Out'); ticklength;
cmap = goodcolormap('bwr',10000); 
colormap(hot); set(gca,'Color','none');

%% Information gain: 1-up 1-down, 10-up 10-down
load('staircases_vs_stepsize.mat');
subplot(3,3,3); hold on;
M = 1; N = 1;
x_stim1 = squeeze(ConsecutiveStaircase(2).StationaryDistribution.x_stim(M,N,:));
y_stim1 = squeeze(ConsecutiveStaircase(2).StationaryDistribution.y(M,N,:));
plot(x_stim1,y_stim1,'--k');
x_stim2 = squeeze(ConsecutiveStaircase(4).StationaryDistribution.x_stim(M,N,:));
y_stim2 = interp1(x_stim2,squeeze(ConsecutiveStaircase(4).StationaryDistribution.y(M,N,:)),x_stim1);
y_stim2 = y_stim2/sum(y_stim2); % re-normalise
plot(x_stim1,y_stim2,'k');
ylim([0 0.15]); set(gca,'XTick',[-2 0 2],'XTickLabel',{'-2\sigma','0','+2\sigma'},'YTick',0:0.05:0.15);
xlabel('Stimulus, \it s'); ylabel('Stationary probability, \pi(s)');
legend('0.05\sigma','0.2\sigma','Color','none'); hleg = legend('show'); title(hleg,'Step-size');
set(gca,'TickDir','Out'); ticklength;
set(gca,'Color','none');

%% Information gain vs Step-size
subplot(3,3,6); hold on;
M = 1; N = 1;
for k=1:numel(NonconsecutiveStaircase)
    StepSize(k) = NonconsecutiveStaircase(k).DesignParameter.stepsize;
    InfoGain(k) = NonconsecutiveStaircase(k).StationaryDistribution.EntropyUniform(M,N) - NonconsecutiveStaircase(k).StationaryDistribution.Entropy(M,N);
end
plot(StepSize,InfoGain,'-k'); axis([0.005 1 0 5.35]); set(gca,'XScale','Log','XTick',[0.01 0.1 1],'XTickLabel',{'0.01\sigma','0.1\sigma','1\sigma'}); 
set(gca,'TickDir','Out'); ticklength;
xlabel('Step-size'); ylabel('Information gain [bits]');
set(gca,'Color','none');