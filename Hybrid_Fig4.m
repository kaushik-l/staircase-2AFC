%% Compute asymptotic accuracy for different hybrid staircases
MNmax = 2; gamma = 0.01:0.01:5;
for i=1:length(gamma)
    accuracy_cont(:,:,i) = 100*ComputeAccuracyContinuous(MNmax,gamma(i));
    accuracy_blkd(:,:,i) = 100*ComputeAccuracyBlocked(MNmax,gamma(i));
    accuracy_noncons(:,:,i) = 100*ComputeAccuracyNonConsec(MNmax,gamma(i));
end

figure; hold on; set(gcf,'Position',[1107 1037 900 900]);

subplot(2,2,1); hold on;
%% Weighted Continuous
plot(gamma,squeeze(accuracy_cont(1,2,:)),'Color',[0 0.25 1],'linewidth',2);

%% Weighted Blocked
plot(gamma,squeeze(accuracy_blkd(1,2,:)),'r','linewidth',2);

%% Weighted Non-consecutive
plot(gamma,squeeze(accuracy_noncons(1,2,:)),'Color',[1 0.5 0],'linewidth',2);

%%
axis([0 5 50 100]); set(gca,'XTick',0:1:5,'YTick',50:10:100);
legend({'Weighted Continuous','Weighted Blocked','Weighted Non-consecutive'},'Location','nw');
xlabel('Weight, \it {\gamma = \Delta_u /\Delta_d}','fontsize',15);
ylabel('Asymptotic Accuracy','fontsize',15);