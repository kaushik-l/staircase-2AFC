% load staircase solutions
load('staircases.mat');

%% define colors
clrs = [0 0.4 0.8 ; 0 0.8 0.4 ; 0.8 0 0.8 ; 0.8 0.4 0];
figure; hold on; set(gcf,'Position',[50 300 1000 1000]);

%% 1-up 1-down, 10-up 10-down
subplot(3,3,1); hold on;
M = 1; N = 3; stim_max = ConsecutiveStaircase.DesignParameter.stim_max;
imagesc(squeeze(ConsecutiveStaircase.ConvergenceRate.x_stim(M,N,:)),...
    squeeze(ConsecutiveStaircase.ConvergenceRate.x_stim(M,N,:)),...
    squeeze(ConsecutiveStaircase.ConvergenceRate.TPM(M,N,:,:)),[0 1]);
axis([-stim_max stim_max -stim_max stim_max]); 
set(gca,'XTick',[-2 0 2],'XTickLabel',{'-2\sigma','0','+2\sigma'},'YTick',[-2 0 2],'YTickLabel',{'-2\sigma','0','+2\sigma'});
set(gca,'TickDir','Out'); ticklength;
cmap = goodcolormap('wr',10000); colormap(gca,cmap); 
cb = colorbar('northoutside','AxisLocation','in'); cb.Position = cb.Position + [0.05 0.07 -0.1 0]; 
cb.Label.String = 'Transition Prob.'; cb.Label.FontSize = 10; cb.Label.Position = cb.Label.Position + [0 4 0];
title('1-up N-down','fontweight','bold'); xlabel('Stimulus, \it s'); ylabel('Stimulus, \it s'); 

%% Error vs trials for different (M,N)
subplot(3,3,2); hold on;
maxtrials = 1e4;
M = 1; N = 1;
plot(1:maxtrials,ConsecutiveStaircase.ConvergenceRate.lambda(M,N).^(1:maxtrials),'Color',clrs(2,:));
M = 10; N = 10;
plot(1:maxtrials,ConsecutiveStaircase.ConvergenceRate.lambda(M,N).^(1:maxtrials),'--','Color',clrs(2,:));
set(gca,'XScale','Log','YTick',0:0.25:1); xlim([1 maxtrials]);
title('Convergence rate','fontweight','bold'); xlabel('Number of trials, \it t'); ylabel('Distance from balance point, {\epsilon_{\itt} /\epsilon_0}'); 

%% Number of trials vs N
subplot(3,3,3); hold on;
plot(1:NonconsecutiveStaircase.DesignParameter.N_max,diag(NonconsecutiveStaircase.ConvergenceRate.numtrials),'-','Color',clrs(1,:));
plot(1:ConsecutiveStaircase.DesignParameter.N_max,diag(ConsecutiveStaircase.ConvergenceRate.numtrials),'-','Color',clrs(2,:));
plot(1:BlockedStaircase.DesignParameter.N_max,diag(BlockedStaircase.ConvergenceRate.numtrials),'-','Color',clrs(3,:));
plot(1:ContinuousStaircase.DesignParameter.N_max,diag(ContinuousStaircase.ConvergenceRate.numtrials),'-','Color',clrs(4,:));
legend('Nonconsecutive', 'Consecutive', 'Blocked', 'Continuous', 'Location','nw'); hleg = legend('show'); title(hleg,'Update rule');
set(gca,'YScale','Log'); axis([1 10 10 1e5]);
set(gca,'TickDir','Out'); ticklength;
title('99% convergence','fontweight','bold'); xlabel('\it N'); ylabel('Trials to convergence, \it t*');

%% Number of trials vs (M,N)
clr_range = [1 5];
subplot(3,3,4); hold on;
imagesc(1:NonconsecutiveStaircase.DesignParameter.M_max,1:NonconsecutiveStaircase.DesignParameter.N_max,...
    log10(NonconsecutiveStaircase.ConvergenceRate.numtrials),clr_range);
axis([0.5 10.5 0.5 10.5]); set(gca,'YDir','normal'); xlabel('\it M'); ylabel({'\it N'}); title('Nonconsecutive','fontweight','bold');
set(gca,'TickDir','Out'); ticklength;
colormap(gca,hot);
subplot(3,3,5); hold on;
imagesc(1:ConsecutiveStaircase.DesignParameter.M_max,1:ConsecutiveStaircase.DesignParameter.N_max,...
    log10(ConsecutiveStaircase.ConvergenceRate.numtrials),clr_range);
axis([0.5 10.5 0.5 10.5]); set(gca,'YDir','normal'); xlabel('\it M'); ylabel({'\it N'}); title('Consecutive','fontweight','bold');
set(gca,'TickDir','Out'); ticklength;
colormap(gca,hot);
subplot(3,3,7); hold on;
imagesc(1:BlockedStaircase.DesignParameter.M_max,1:BlockedStaircase.DesignParameter.N_max,...
    log10(BlockedStaircase.ConvergenceRate.numtrials),clr_range);
axis([0.5 10.5 0.5 10.5]); set(gca,'YDir','normal'); xlabel('\it M'); ylabel({'\it N'}); title('Blocked','fontweight','bold');
set(gca,'TickDir','Out'); ticklength;
colormap(gca,hot);
cb = colorbar('northoutside','AxisLocation','in','Ticks',[1 3 5]); cb.Position = cb.Position + [0.2 0.075 -0.1 0]; 
cb.Label.String = 'Log_{10}\it t*'; cb.Label.FontSize = 10; cb.Label.Position = cb.Label.Position + [0 4 0];
subplot(3,3,8); hold on;
imagesc(1:ContinuousStaircase.DesignParameter.M_max,1:ContinuousStaircase.DesignParameter.N_max,...
    log10(ContinuousStaircase.ConvergenceRate.numtrials),clr_range);
axis([0.5 10.5 0.5 10.5]); set(gca,'YDir','normal'); xlabel('\it M'); ylabel({'\it N'}); title('Continuous','fontweight','bold');
set(gca,'TickDir','Out'); ticklength;
colormap(gca,hot);

%% Number of trials: 1-up 1-down, 10-up 10-down
load('staircases_vs_stepsize.mat');
subplot(3,3,6); hold on;
maxtrials = 500;
M = 1; N = 1;
plot(1:maxtrials,NonconsecutiveStaircase(2).ConvergenceRate.lambda(M,N).^(1:maxtrials),'--','Color',clrs(2,:));
plot(1:maxtrials,NonconsecutiveStaircase(4).ConvergenceRate.lambda(M,N).^(1:maxtrials),'Color',clrs(2,:));
set(gca,'XScale','Log','YTick',0:0.25:1); xlim([1 maxtrials]);
set(gca,'TickDir','Out'); ticklength;
xlabel('Number of trials, \it t'); ylabel('Distance from balance point, {\epsilon_{\itt} /\epsilon_0}'); 
legend('0.05\sigma','0.2\sigma'); hleg = legend('show'); title(hleg,'Step-size');

%% Number of trials vs Step-size
subplot(3,3,9); hold on;
M = 1; N = 1;
for k=1:numel(NonconsecutiveStaircase)
    StepSize(k) = NonconsecutiveStaircase(k).DesignParameter.stepsize;
    NumTrials(k) = NonconsecutiveStaircase(k).ConvergenceRate.numtrials(M,N);
end
plot(StepSize,NumTrials,'-k'); axis([0.005 1 1 1e3]);
set(gca,'XScale','Log','XTick',[0.01 0.1 1],'XTickLabel',{'0.01\sigma','0.1\sigma','1\sigma'},'YScale','Log');
set(gca,'TickDir','Out'); ticklength;
xlabel('Step-size'); ylabel('Trials to convergence, \it t*');