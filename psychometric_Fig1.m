%% psychometric function
x = linspace(-5,5,1000);
f = normpdf(x,0,2); F = normcdf(x,0,2);
figure; hold on; title('Figure 1');
subplot(2,4,1); plot(x,F); vline(0,'-k'); box off; %axis off;
xlabel('Stimulus, \it s','Interpreter','tex'); ylabel({'Percentage of'; 'rightward response'});
F = 0.5 + abs(F-0.5);
subplot(2,4,5); hold on; plot(x,F); axis([0 5 0.5 1]); vline(0,'-k'); box off; %axis off;
xlabel({'Stimulus deviation, ' ; '\deltas = |\it s - \it s_0|'},'Interpreter','tex'); ylabel({'Percentage of'; 'correct response'});

%% 1-up 1-down probs
N = 1; M = 1;
Pd = F.^N; Pu = (1-F).^M;
subplot(2,4,2); hold on; plot(x,F); axis([0 5 0.5 1]); vline(0,'-k'); box off; ylabel({'Percentage of'; 'rightward response'}); %axis off;
indx = find((abs(Pu - Pd)) == min(abs(Pu - Pd)),1,'last'); plot(x(indx),F(indx),'ok','MarkerFaceColor','k');
subplot(2,4,6); hold on; plot(x,Pd); plot(x,Pu); axis([0 5 0 1]); vline(0,'-k'); box off; %axis off;
title('1-up 1-down','fontweight','bold');

%% 1-up 2-down probs
N = 2; M = 1;
Pd = F.^N; Pu = (1-F).^M;
subplot(2,4,3); hold on; plot(x,F); axis([0 5 0.5 1]); vline(0,'-k'); box off; %axis off;
indx = find((abs(Pu - Pd)) == min(abs(Pu - Pd)),1,'last'); plot(x(indx),F(indx),'ok','MarkerFaceColor','k');
subplot(2,4,7); hold on; plot(x,Pd); plot(x,Pu); axis([0 5 0 1]); vline(0,'-k'); box off; %axis off;
title('1-up 2-down','fontweight','bold');

%% 1-up 10-down probs
N = 10; M = 1;
Pd = F.^N; Pu = (1-F).^M;
subplot(2,4,4); hold on; plot(x,F); axis([0 5 0.5 1]); vline(0,'-k'); box off; %axis off;
indx = find((abs(Pu - Pd)) == min(abs(Pu - Pd)),1,'last'); plot(x(indx),F(indx),'ok','MarkerFaceColor','k');
subplot(2,4,8); hold on; plot(x,Pd); plot(x,Pu); axis([0 5 0 1]); vline(0,'-k'); box off; %axis off;
title('1-up 10-down','fontweight','bold');

%%
sgtitle('Figure 1');