%% Compute asymptotic accuracy for different staircases
MNmax = 10;
accuracy_cont = 100*ComputeAccuracyContinuous(MNmax);
accuracy_blkd = 100*ComputeAccuracyBlocked(MNmax);
accuracy_noncons = 100*ComputeAccuracyNonConsec(MNmax);

figure; hold on; set(gcf,'Position',[1107 1037 1690 934]); title('Figure 2');

%% 1-up N-down
subplot(2,3,1); hold on;
plot(1:MNmax,accuracy_cont(1,:),'Color',[0 0.25 1]);
plot(1:MNmax,accuracy_cont(1,:),'o','Color',[0 0.25 1],'MarkerFaceColor',[0 0.25 1]);
plot(1:MNmax,accuracy_blkd(1,:),'Color','r');
plot(1:MNmax,accuracy_blkd(1,:),'o','Color','r','MarkerFaceColor','r');
plot(1:MNmax,accuracy_noncons(1,:),'Color',[1 0.5 0]);
plot(1:MNmax,accuracy_noncons(1,:),'o','Color',[1 0.5 0],'MarkerFaceColor',[1 0.5 0]);
axis([1 10 50 100]);
title('1-up N-down','fontsize',20);

%% M-up 10-down
subplot(2,3,2); hold on;
plot(1:MNmax,accuracy_cont(:,MNmax),'Color',[0 0.25 1]);
plot(1:MNmax,accuracy_cont(:,MNmax),'o','Color',[0 0.25 1],'MarkerFaceColor',[0 0.25 1]);
plot(1:MNmax,accuracy_blkd(:,MNmax),'Color','r');
plot(1:MNmax,accuracy_blkd(:,MNmax),'o','Color','r','MarkerFaceColor','r');
plot(1:MNmax,accuracy_noncons(:,MNmax),'Color',[1 0.5 0]);
plot(1:MNmax,accuracy_noncons(:,MNmax),'o','Color',[1 0.5 0],'MarkerFaceColor',[1 0.5 0]);
axis([1 10 50 100]);
title('M-up 10-down','fontsize',20);

%% M-up 2M-down
accuracy_cont = 100*ComputeAccuracyContinuous(2*MNmax);
accuracy_blkd = 100*ComputeAccuracyBlocked(2*MNmax);
accuracy_noncons = 100*ComputeAccuracyNonConsec(2*MNmax);
subplot(2,3,3); hold on;
plot(1:MNmax,diag(accuracy_cont(1:MNmax,2*(1:MNmax))),'Color',[0 0.25 1]);
plot(1:MNmax,diag(accuracy_cont(1:MNmax,2*(1:MNmax))),'o','Color',[0 0.25 1],'MarkerFaceColor',[0 0.25 1]);
plot(1:MNmax,diag(accuracy_blkd(1:MNmax,2*(1:MNmax))),'Color','r');
plot(1:MNmax,diag(accuracy_blkd(1:MNmax,2*(1:MNmax))),'o','Color','r','MarkerFaceColor','r');
plot(1:MNmax,diag(accuracy_noncons(1:MNmax,2*(1:MNmax))),'Color',[1 0.5 0]);
plot(1:MNmax,diag(accuracy_noncons(1:MNmax,2*(1:MNmax))),'o','Color',[1 0.5 0],'MarkerFaceColor',[1 0.5 0]);
axis([1 10 60 75]);
title('M-up N-down','fontsize',20);

%%
accuracy_cont = ComputeAccuracyContinuous(MNmax);
accuracy_blkd = ComputeAccuracyBlocked(MNmax);
accuracy_noncons = ComputeAccuracyNonConsec(MNmax);
%% Continuous staircase
subplot(2,3,4); hold on;
% draw circles proportional to accuracy
for m = 1:MNmax
    for n = 1:MNmax
        if m < n
            plot(n,m,'ro','markersize',20);
            plot(n,m,'bo','markersize',20*(2*accuracy_cont(m,n)-1),'markerfacecolor','b');
        end
    end
end
% draw a triangle around the region m > n
hold on;
plot(1:MNmax,1:MNmax,'k');
plot(1:MNmax,MNmax*ones(1,MNmax),'k');
plot(ones(1,MNmax),1:MNmax,'k');
text(1.5,0.7*MNmax,'Chance (50%)','fontsize',16);
% labels
title('Asymptotic Accuracy (Continuous)','fontsize',20);
xlabel('N','fontsize',14), ylabel('M','fontsize',14);
set(gca,'fontsize',14);
axis([0 MNmax+1 0 MNmax+1]);
axis square;

%% Blocked staircase
subplot(2,3,5); hold on;
% draw circles proportional to accuracy
for m = 1:MNmax
    for n = 1:MNmax
        if m < n
            plot(n,m,'ro','markersize',20);
            plot(n,m,'bo','markersize',20*(2*accuracy_blkd(m,n)-1),'markerfacecolor','b');
        end
    end
end
% draw a triangle around the region m > n
hold on;
plot(1:MNmax,1:MNmax,'k');
plot(1:MNmax,MNmax*ones(1,MNmax),'k');
plot(ones(1,MNmax),1:MNmax,'k');
text(1.5,0.7*MNmax,'Chance (50%)','fontsize',16);
% labels
title('Asymptotic Accuracy (Blocked)','fontsize',20);
xlabel('N','fontsize',14), ylabel('M','fontsize',14);
set(gca,'fontsize',14);
axis([0 MNmax+1 0 MNmax+1]);
axis square;

%% Non-consecutive staircase
subplot(2,3,6); hold on;
% draw circles proportional to accuracy
for m = 1:MNmax
    for n = 1:MNmax
        if m < n
            plot(n,m,'ro','markersize',20);
            plot(n,m,'bo','markersize',20*(2*accuracy_noncons(m,n)-1),'markerfacecolor','b');
        end
    end
end
% draw a triangle around the region m > n
hold on;
plot(1:MNmax,1:MNmax,'k');
plot(1:MNmax,MNmax*ones(1,MNmax),'k');
plot(ones(1,MNmax),1:MNmax,'k');
text(1.5,0.7*MNmax,'Chance (50%)','fontsize',16);
% labels
title('Asymptotic Accuracy (Non-consecutive)','fontsize',20);
xlabel('N','fontsize',14), ylabel('M','fontsize',14);
set(gca,'fontsize',14);
axis([0 MNmax+1 0 MNmax+1]);
axis square;

%%
sgtitle('Figure 2');