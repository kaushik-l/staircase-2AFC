figure; hold on; set(gcf,'Position',[1107 1037 1200 800]);

%% 1-up 2-down staircase (effect of update method)
s0 = 3; ds = 0.1; ntrls = 200; M = 1; N = 2; error = 0.01;
[curve_cont,trls_cont,curve_blkd,trls_blkd,curve_noncons,trls_noncons] = SimulateStaircase(s0,ds,ntrls,M,N,error);
subplot(2,3,1); hold on;
plot(1:ntrls+1,curve_cont,'Color',[0 0.25 1],'linewidth',2);
plot(1:ntrls+1,curve_blkd,'Color','r','linewidth',2);
plot(1:ntrls+1,curve_noncons,'Color',[1 0.5 0],'linewidth',2);
axis([0 200 0.1 s0]); set(gca,'YScale','Log','YTick',[0.1 1 2 3],'XTick',50:50:200);
legend({'Continuous','Blocked','Non-consecutive'});
xlabel('Number of trials'); ylabel({'Stimlulus amplitude' ; '(in units of \sigma)'});

%% 1-up N-down (effect of N)
clear;
s0 = 3; ds = 0.1; ntrls = 300; M = 1; Nmax = 10; error = 0.01;
for N=1:Nmax
    [~,trls_cont(N),~,trls_blkd(N),~,trls_noncons(N)] = SimulateStaircase(s0,ds,ntrls,M,N,error);
end
subplot(2,3,2); hold on;
plot(1:Nmax,trls_cont,'Marker','o','Color',[0 0.25 1],'MarkerFaceColor',[0 0.25 1],'linewidth',2);
plot(1:Nmax,trls_blkd,'Marker','o','Color','r','MarkerFaceColor','r','linewidth',2);
plot(1:Nmax,trls_noncons,'Marker','o','Color',[1 0.5 0],'MarkerFaceColor',[1 0.5 0],'linewidth',2);
axis([1 Nmax 50 300]); set(gca,'YScale','Log','YTick',50:50:300,'XTick',2:2:10);
xlabel('\it N'); ylabel('Number of trials');

%% M-up 10-down (effect of M)
clear;
s0 = 3; ds = 0.1; ntrls = 100000; Mmax = 10; N = 10; error = 0.01;
for M=1:Mmax
    [~,trls_cont(M),~,trls_blkd(M),~,trls_noncons(M)] = SimulateStaircase(s0,ds,ntrls,M,N,error);
end
subplot(2,3,3); hold on;
plot(1:Mmax,trls_cont,'Marker','o','Color',[0 0.25 1],'MarkerFaceColor',[0 0.25 1],'linewidth',2);
plot(1:Mmax,trls_blkd,'Marker','o','Color','r','MarkerFaceColor','r','linewidth',2);
plot(1:Mmax,trls_noncons,'Marker','o','Color',[1 0.5 0],'MarkerFaceColor',[1 0.5 0],'linewidth',2);
axis([1 Mmax 1e1 1e4]); set(gca,'YScale','Log','YTick',[1e2 1e3 1e4],'XTick',2:2:10);
xlabel('\it M');

%% 1-up 2-down (effect of ds)
clear;
s0 = 3; dsmax = 0.5; ntrls = 10000; M = 1; N = 2; error = 0.01; count = 0;
ds = [0.01 0.1:0.1:dsmax];
for ds_iter=ds
    count = count+1;
    [~,trls_cont(count),~,trls_blkd(count),~,trls_noncons(count)] = SimulateStaircase(s0,ds_iter,ntrls,M,N,error);
end
subplot(2,3,4); hold on;
plot(ds,1./trls_cont,'Color',[0 0.25 1],'MarkerFaceColor',[0 0.25 1],'linewidth',2);
plot(ds,1./trls_blkd,'Color','r','MarkerFaceColor','r','linewidth',2);
plot(ds,1./trls_noncons,'Color',[1 0.5 0],'MarkerFaceColor',[1 0.5 0],'linewidth',2);
axis([0 0.5 0 0.1]); set(gca,'YTick',0:0.02:0.1,'XTick',0:0.1:0.5);
xlabel('Step size, {\it \Delta} (in units of \sigma)'); ylabel({'Convergence rate'; '(No. of trials)^{-1}'});

%% 1-up 2-down (effect of s0)
clear;
MNmax = 2;
s0max = 2; ds = 0.1; ntrls = 200; M = 1; N = 2; error = 0.01;
% continuous
accuracy_cont = ComputeAccuracyContinuous(MNmax); accuracy_cont = accuracy_cont(1,2);
s0 = accuracy_cont + [0:0.1:s0max]; count = 0;
for s0_iter=s0
    count = count+1;
    [~,trls_cont(count),~,~,~,~] = SimulateStaircase(s0_iter,ds,ntrls,M,N,error);
end
% blocked
accuracy_blkd = ComputeAccuracyBlocked(MNmax); accuracy_blkd = accuracy_blkd(1,2);
s0 = accuracy_blkd + [0:0.1:s0max]; count = 0;
for s0_iter=s0
    count = count+1;
    [~,~,~,trls_blkd(count),~,~] = SimulateStaircase(s0_iter,ds,ntrls,M,N,error);
end
% non-consecutive
accuracy_noncons = ComputeAccuracyNonConsec(MNmax); accuracy_noncons = accuracy_noncons(1,2);
s0 = accuracy_noncons + [0:0.1:s0max]; count = 0;
for s0_iter=s0
    count = count+1;
    [~,~,~,~,~,trls_noncons(count)] = SimulateStaircase(s0_iter,ds,ntrls,M,N,error);
end
subplot(2,3,5); hold on;
plot(0:0.1:s0max,trls_cont,'Color',[0 0.25 1],'MarkerFaceColor',[0 0.25 1],'linewidth',2);
plot(0:0.1:s0max,trls_blkd,'Color','r','MarkerFaceColor','r','linewidth',2);
plot(0:0.1:s0max,trls_noncons,'Color',[1 0.5 0],'MarkerFaceColor',[1 0.5 0],'linewidth',2);
axis([0 0.5 0 0.1]); set(gca,'YTick',0:0.02:0.1,'XTick',0:0.1:0.5);
% xlabel('Step size, {\it \Delta} (in units of \sigma)'); ylabel({'Convergence rate'; '(No. of trials)^{-1}'});