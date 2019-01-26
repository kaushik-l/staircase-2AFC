%% psychometric function
x = linspace(-5,5,1000);
f = normpdf(x,0,2); F = normcdf(x,0,2);
figure; plot(x,F); vline(0,'-k'); box off; %axis off;
F = 0.5 + abs(F-0.5);
figure; plot(x,F); axis([0 5 0.5 1]); vline(0,'-k'); box off; %axis off;

%% 1-up 1-down probs
N = 1; M = 1;
Pd = F.^N; Pu = (1-F).^M;
figure; hold on; plot(x,Pd); plot(x,Pu); axis([0 5 0 1]); vline(0,'-k'); box off; %axis off;

%% 1-up 2-down probs
N = 2; M = 1;
Pd = F.^N; Pu = (1-F).^M;
figure; hold on; plot(x,Pd); plot(x,Pu); axis([0 5 0 1]); vline(0,'-k'); box off; %axis off;

%% 1-up 10-down probs
N = 10; M = 1;
Pd = F.^N; Pu = (1-F).^M;
figure; hold on; plot(x,Pd); plot(x,Pu); axis([0 5 0 1]); vline(0,'-k'); box off; %axis off;