% load staircase solutions
load('staircases.mat');

%% Psychometric function
x = linspace(-5,5,1000); sigma = 2;
f = normpdf(x,0,sigma); F = normcdf(x,0,sigma);
figure; hold on; set(gcf,'Position',[50 300 1000 600]);

%% 1-up 1-down probs
% Draw Pup and Pdown
M = 1; N = 1; 
Pd = F/N; Pu = (1-F)/M;
subplot(2,3,1); hold on; plot(x,Pd,'LineWidth',2); plot(x,Pu,'LineWidth',2); vline(sigma*norminv(NonconsecutiveStaircase.BalancePoint.Percentile(M,N)),'--k'); 
box off; set(gca,'YTick',0.25:0.25:1,'XTick',sigma*[-2.5 0 2.5],'XTickLabel',{'-2.5\sigma','0','+2.5\sigma'},'FontSize',10); %axis off;
set(gca,'TickDir','Out'); ticklength; set(gca,'Color','none');
title('1-up 1-down','fontweight','bold'); ylabel({'Transition probability'}); 
legend('Downward', 'Upward','Location','east','Color','None'); hleg = legend('show'); title(hleg,'Direction');
% Draw psychometric function
subplot(2,3,4); plot(x,F,'k','LineWidth',2); vline(sigma*norminv(NonconsecutiveStaircase.BalancePoint.Percentile(M,N)),'--k'); 
hline(NonconsecutiveStaircase.BalancePoint.Percentile(M,N),'--k'); box off; %axis off;
xlabel('Stimulus, \it s','Interpreter','tex'); ylabel({'Cumulative Probability, \Phi(s)'}); 
set(gca,'YTick',0.25:0.25:1,'XTick',sigma*[-2.5 0 2.5],'XTickLabel',{'-2.5\sigma','0','+2.5\sigma'},'FontSize',10); %axis off;
set(gca,'TickDir','Out'); ticklength; set(gca,'Color','none');
F_targ = NonconsecutiveStaircase.BalancePoint.Percentile(M,N);
[~,indx] = min(abs(F - F_targ));
hold on; plot(x(indx),F(indx),'ok','MarkerFaceColor','k');

%% 1-up 3-down probs
% Draw Pup and Pdown
M = 1; N = 3; 
Pd = F/N; Pu = (1-F)/M;
subplot(2,3,2); hold on; plot(x,Pd,'LineWidth',2); plot(x,Pu,'LineWidth',2); vline(sigma*norminv(NonconsecutiveStaircase.BalancePoint.Percentile(M,N)),'--k'); 
box off; set(gca,'YTick',0.25:0.25:1,'XTick',sigma*[-2.5 0 2.5],'XTickLabel',{'-2.5\sigma','0','+2.5\sigma'},'FontSize',10); %axis off;
set(gca,'TickDir','Out'); ticklength; set(gca,'Color','none');
title('1-up 3-down','fontweight','bold');
% Draw psychometric function
subplot(2,3,5); plot(x,F,'k','LineWidth',2); vline(sigma*norminv(NonconsecutiveStaircase.BalancePoint.Percentile(M,N)),'--k'); 
hline(NonconsecutiveStaircase.BalancePoint.Percentile(M,N),'--k'); box off; %axis off; box off; %axis off;
xlabel('Stimulus, \it s','Interpreter','tex'); 
set(gca,'YTick',0.25:0.25:1,'XTick',sigma*[-2.5 0 2.5],'XTickLabel',{'-2.5\sigma','0','+2.5\sigma'},'FontSize',10); %axis off;
set(gca,'TickDir','Out'); ticklength; set(gca,'Color','none');
F_targ = NonconsecutiveStaircase.BalancePoint.Percentile(M,N);
[~,indx] = min(abs(F - F_targ));
hold on; plot(x(indx),F(indx),'ok','MarkerFaceColor','k');

%% 3-up 1-down probs
% Draw Pup and Pdown
M = 3; N = 1; 
Pd = F/N; Pu = (1-F)/M;
subplot(2,3,3); hold on; plot(x,Pd,'LineWidth',2); plot(x,Pu,'LineWidth',2); vline(sigma*norminv(NonconsecutiveStaircase.BalancePoint.Percentile(M,N)),'--k'); 
box off; set(gca,'YTick',0.25:0.25:1,'XTick',sigma*[-2.5 0 2.5],'XTickLabel',{'-2.5\sigma','0','+2.5\sigma'},'FontSize',10); %axis off;
set(gca,'TickDir','Out'); ticklength; set(gca,'Color','none');
title('3-up 1-down','fontweight','bold');
% Draw psychometric function
subplot(2,3,6); plot(x,F,'k','LineWidth',2); vline(sigma*norminv(NonconsecutiveStaircase.BalancePoint.Percentile(M,N)),'--k'); 
hline(NonconsecutiveStaircase.BalancePoint.Percentile(M,N),'--k'); box off; %axis off; box off; %axis off;
xlabel('Stimulus, \it s','Interpreter','tex'); 
set(gca,'YTick',0.25:0.25:1,'XTick',sigma*[-2.5 0 2.5],'XTickLabel',{'-2.5\sigma','0','+2.5\sigma'},'FontSize',10); %axis off;
set(gca,'TickDir','Out'); ticklength; set(gca,'Color','none');
F_targ = NonconsecutiveStaircase.BalancePoint.Percentile(M,N);
[~,indx] = min(abs(F - F_targ));
hold on; plot(x(indx),F(indx),'ok','MarkerFaceColor','k');