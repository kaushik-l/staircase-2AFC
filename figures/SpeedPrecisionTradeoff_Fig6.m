clear;
load('staircases_vs_stepsize.mat');

%% define colors
clrs = [0 0.4 0.8 ; 0 0.8 0.4 ; 0.8 0 0.8 ; 0.8 0.4 0];
figure; hold on; set(gcf,'Position',[50 300 600 500]);

%% Information gain vs Step-size
subplot(2,2,1); hold on;
M = 1; N = 1;
for k=1:numel(NonconsecutiveStaircase)
    StepSize(k) = NonconsecutiveStaircase(k).DesignParameter.stepsize;
    InfoGain(k) = NonconsecutiveStaircase(k).StationaryDistribution.EntropyUniform(M,N) - NonconsecutiveStaircase(k).StationaryDistribution.Entropy(M,N);
    NumTrials(k) = NonconsecutiveStaircase(k).ConvergenceRate.numtrials(M,N);
end
yyaxis left;
plot(StepSize,InfoGain); axis([0.01 1 0 5.35]);
ylabel('Information gain [bits]');
yyaxis right;
plot(StepSize,1./NumTrials); axis([0.01 1 1e-3 1]); set(gca,'YScale','Log');
set(gca,'XScale','Log','XTick',[0.01 0.1 1],'XTickLabel',{'0.01\sigma','0.1\sigma','1\sigma'});
set(gca,'TickDir','Out'); ticklength;
xlabel('Step-size');
set(gca,'Color','none');

% subplot(2,2,2); hold on;
% M = 1; N = 1;
% plot(StepSize,InfoGain./NumTrials); %axis([0.01 1 1e-3 1]); set(gca,'YScale','Log');
% ylabel('Rate of Information gain [bits/trial]');
% set(gca,'XScale','Log','XTick',[0.01 0.1 1],'XTickLabel',{'0.01\sigma','0.1\sigma','1\sigma'});
% set(gca,'TickDir','Out'); ticklength;
% xlabel('Step-size');

%% load staircase solutions
load('staircases.mat');

subplot(2,2,2); hold on;
yyaxis left;
plot(1:ConsecutiveStaircase.DesignParameter.N_max,...
    diag(ConsecutiveStaircase.StationaryDistribution.EntropyUniform) - diag(ConsecutiveStaircase.StationaryDistribution.Entropy));
yyaxis right;
plot(1:ConsecutiveStaircase.DesignParameter.N_max,...
    1./diag(ConsecutiveStaircase.ConvergenceRate.numtrials)); set(gca,'YScale','Log');
xlim([1 ConsecutiveStaircase.DesignParameter.N_max]);
xlabel('\it N,M'); ylabel('Convergence rate, \it{^1/_{t*}} \rm [trial^{-1}]');
set(gca,'TickDir','Out'); ticklength;
set(gca,'Color','none');

%% 
subplot(2,2,3); hold on;
M = 1;
NumTrials = NonconsecutiveStaircase.ConvergenceRate.numtrials(M,:);
InfoGain = NonconsecutiveStaircase.StationaryDistribution.EntropyUniform(M,:) - NonconsecutiveStaircase.StationaryDistribution.Entropy(M,:);
plot(NonconsecutiveStaircase.BalancePoint.Percentile(M,:), InfoGain./NumTrials, 'Color', clrs(1,:));

NumTrials = ConsecutiveStaircase.ConvergenceRate.numtrials(M,:);
InfoGain = ConsecutiveStaircase.StationaryDistribution.EntropyUniform(M,:) - ConsecutiveStaircase.StationaryDistribution.Entropy(M,:);
plot(ConsecutiveStaircase.BalancePoint.Percentile(M,:), InfoGain./NumTrials, 'Color', clrs(2,:));

NumTrials = BlockedStaircase.ConvergenceRate.numtrials(M,:);
InfoGain = BlockedStaircase.StationaryDistribution.EntropyUniform(M,:) - BlockedStaircase.StationaryDistribution.Entropy(M,:);
plot(BlockedStaircase.BalancePoint.Percentile(M,:), InfoGain./NumTrials, 'Color', clrs(3,:));

NumTrials = ContinuousStaircase.ConvergenceRate.numtrials(M,:);
InfoGain = ContinuousStaircase.StationaryDistribution.EntropyUniform(M,:) - ContinuousStaircase.StationaryDistribution.Entropy(M,:);
plot(ContinuousStaircase.BalancePoint.Percentile(M,:), InfoGain./NumTrials, 'Color', clrs(4,:));

N = 1;
NumTrials = NonconsecutiveStaircase.ConvergenceRate.numtrials(:,N);
InfoGain = NonconsecutiveStaircase.StationaryDistribution.EntropyUniform(:,N) - NonconsecutiveStaircase.StationaryDistribution.Entropy(:,N);
plot(NonconsecutiveStaircase.BalancePoint.Percentile(:,N), InfoGain./NumTrials, 'Color', clrs(1,:));

NumTrials = ConsecutiveStaircase.ConvergenceRate.numtrials(:,N);
InfoGain = ConsecutiveStaircase.StationaryDistribution.EntropyUniform(:,N) - ConsecutiveStaircase.StationaryDistribution.Entropy(:,N);
plot(ConsecutiveStaircase.BalancePoint.Percentile(:,N), InfoGain./NumTrials, 'Color', clrs(2,:));

NumTrials = BlockedStaircase.ConvergenceRate.numtrials(:,N);
InfoGain = BlockedStaircase.StationaryDistribution.EntropyUniform(:,N) - BlockedStaircase.StationaryDistribution.Entropy(:,N);
plot(BlockedStaircase.BalancePoint.Percentile(:,N), InfoGain./NumTrials, 'Color', clrs(3,:));

NumTrials = ContinuousStaircase.ConvergenceRate.numtrials(:,N);
InfoGain = ContinuousStaircase.StationaryDistribution.EntropyUniform(:,N) - ContinuousStaircase.StationaryDistribution.Entropy(:,N);
plot(ContinuousStaircase.BalancePoint.Percentile(:,N), InfoGain./NumTrials, 'Color', clrs(4,:));

xlabel({'Cumulative prob. of' ; 'balance point, \Phi(s*)'}); ylabel('Efficiency [bits/trial]');
set(gca,'YTick',0:0.01:0.05);
set(gca,'TickDir','Out'); ticklength;
legend('Nonconsecutive', 'Consecutive', 'Blocked', 'Continuous', 'Location','nw','Color','None'); hleg = legend('show'); title(hleg,'Update algorithm');
set(gca,'Color','none');

%%
subplot(2,2,4); hold on;
NumTrials = NonconsecutiveStaircase.ConvergenceRate.numtrials;
InfoGain = NonconsecutiveStaircase.StationaryDistribution.EntropyUniform - NonconsecutiveStaircase.StationaryDistribution.Entropy;
Efficiency(1,:) = InfoGain(:)./NumTrials(:);

NumTrials = ConsecutiveStaircase.ConvergenceRate.numtrials;
InfoGain = ConsecutiveStaircase.StationaryDistribution.EntropyUniform - ConsecutiveStaircase.StationaryDistribution.Entropy;
Efficiency(2,:) = InfoGain(:)./NumTrials(:);

NumTrials = BlockedStaircase.ConvergenceRate.numtrials;
InfoGain = BlockedStaircase.StationaryDistribution.EntropyUniform - BlockedStaircase.StationaryDistribution.Entropy;
Efficiency(3,:) = InfoGain(:)./NumTrials(:);

NumTrials = ContinuousStaircase.ConvergenceRate.numtrials;
InfoGain = ContinuousStaircase.StationaryDistribution.EntropyUniform - ContinuousStaircase.StationaryDistribution.Entropy;
Efficiency(4,:) = InfoGain(:)./NumTrials(:);

MeanEfficiency = mean(Efficiency,2);
SEMEfficiency = std(Efficiency,[],2)/sqrt(100);

for k=1:4, h = bar(k,MeanEfficiency(k),'FaceColor',clrs(k,:)); end
er = errorbar(1:4,MeanEfficiency,SEMEfficiency,'.k'); er.CapSize = 0; 
ylim([0 0.03]); set(gca,'XTick',[],'YTick',0:0.01:0.03); xlabel('Update algorithm');
set(gca,'Color','none');