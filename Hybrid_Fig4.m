%% Compute asymptotic accuracy for different hybrid staircases
MNmax = 10; gamma = 0.01:0.01:2;
for i=1:length(gamma)
    accuracy_cont(:,:,i) = 100*ComputeAccuracyContinuous(MNmax,gamma(i));
    accuracy_blkd(:,:,i) = 100*ComputeAccuracyBlocked(MNmax,gamma(i));
    accuracy_noncons(:,:,i) = 100*ComputeAccuracyNonConsec(MNmax,gamma(i));
end

figure; hold on; set(gcf,'Position',[-1807 300 900 900]);

subplot(2,2,1); hold on;
%% Weighted Continuous
plot(gamma,squeeze(accuracy_cont(1,2,:)),'Color',[0 0.25 1],'linewidth',2);
%% Weighted Blocked
plot(gamma,squeeze(accuracy_blkd(1,2,:)),'r','linewidth',2);
%% Weighted Non-consecutive
plot(gamma,squeeze(accuracy_noncons(1,2,:)),'Color',[1 0.5 0],'linewidth',2);
%% legend
axis([0 2 50 100]); set(gca,'XTick',0:1:2,'YTick',50:10:100);
legend({'Weighted Continuous','Weighted Blocked','Weighted Non-consecutive'},'Location','nw');
xlabel('Weight, \it {\gamma = \Delta_u /\Delta_d}','fontsize',15);
ylabel('Asymptotic Accuracy','fontsize',15);

%% Weighted Continuous
subplot(2,2,2); hold on;
cmap = [0 0 0; 0.25 0.25 0.25; 0.5 0.5 0.5; 0.75 0.75 0.75];
levels = [55 65 75 85];
[C,h] = contour(gamma,1:10,squeeze(accuracy_cont(1,:,:)),levels,'linewidth',2); set(gca,'YDir','Normal');
colormap(cmap);
clabel(C,h,levels,'Fontsize',10,'Color','r'); set(gca,'XTick',0:0.5:2,'YTick',2:2:10);
xlabel('\it \gamma'); ylabel('\it N/M');

%% Weighted Blocked
subplot(2,2,3); hold on;
cmap = [0 0 0; 0.25 0.25 0.25; 0.5 0.5 0.5; 0.75 0.75 0.75];
levels = [55 65 75 85];
[C,h] = contour(gamma,1:10,squeeze(accuracy_blkd(1,:,:)),levels,'linewidth',2); set(gca,'YDir','Normal');
colormap(cmap);
clabel(C,h,levels,'Fontsize',10,'Color','r'); set(gca,'XTick',0:0.5:2,'YTick',2:2:10);
xlabel('\it \gamma'); ylabel('\it N/M');

%% Weighted Non-consecutive
subplot(2,2,4); hold on;
cmap = [0 0 0; 0.25 0.25 0.25; 0.5 0.5 0.5; 0.75 0.75 0.75];
levels = [55 65 75 85];
[C,h] = contour(gamma,1:10,squeeze(accuracy_noncons(1,:,:)),levels,'linewidth',2); set(gca,'YDir','Normal');
colormap(cmap);
clabel(C,h,levels,'Fontsize',10,'Color','r'); set(gca,'XTick',0:0.5:2,'YTick',2:2:10);
xlabel('\it \gamma'); ylabel('\it N/M');
