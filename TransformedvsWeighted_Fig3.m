%% Compute asymptotic accuracy for different transformed and weighted staircase
MNmax = 10; gammamax = 10;
accuracy_cont = 100*ComputeAccuracyContinuous(MNmax);
accuracy_noncons = 100*ComputeAccuracyNonConsec(MNmax);
accuracy_wtd = 100*ComputeAccuracyWeighted(gammamax);

figure; hold on; set(gcf,'Position',[-1807 300 600 600]);

%% Weighted up-down method
p1 = plot(1:gammamax,accuracy_wtd,'Color','k'); axis([1 10 50 100]); set(gca,'XTick',2:2:10); set(gca,'YTick',50:10:100);
xlabel('Weighted up-down, \it {\Delta_u /\Delta_d = \gamma}','fontsize',15);

%% Continuous up-down method
ax1 = gca; ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','left','Color','none'); hold on;
p2 = plot(1:MNmax,accuracy_cont(1,:),'Parent',ax2,'Marker','o','Color',[0 0.25 1],'MarkerFaceColor',[0 0.25 1],'linestyle','none'); hold on;

%% Non-consecutive up-down method
p3 = plot(1:MNmax,accuracy_noncons(1,:),'Parent',ax2,'Marker','o','Color',[1 0.5 0],'MarkerFaceColor',[1 0.5 0],'linestyle','none');
axis([1 10 50 100]); set(gca,'XTick',2:2:10); set(gca,'YTick',50:10:100);
xlabel('Transformed up-down, \it {N/M}','fontsize',15); ylabel('Asymptotic Accuracy','fontsize',15);

%%
legend([p1 p2 p3], {'Weighted','Continuous','Non-consecutive'});