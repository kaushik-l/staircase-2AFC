
function [trls_cont,trls_block,trls_nc,trls_weighted] = learning_curve_saturation(s0,ds,ntrls,M,N,error)

trls_cont = calcStairThresh(s0,ds,ntrls,M,N,error);
trls_block = BlockStairSim(s0,ds,ntrls,M,N,error);
trls_nc = NonConsecStairSim(s0,ds,ntrls,M,N,error);

gamma = N;
trls_weighted = WeightedStairSim(s0,ds,ntrls,gamma,error);

figure, hold on
plot(1:Nmax,trls_cont)
plot(1:Nmax,trls_block)
plot(1:Nmax,trls_nc)
set(gca,'fontsize',14)
title('Number of Trials To Reach Threshold','fontsize',20)
xlabel('N','fontsize',14)
ylabel('Number of Trials','fontsize',14)
ylim([0 250])
legend('Continuous','Blocked','Non-Consecutive','location','northwest')

figure, hold on
plot(gamma,trls_weighted)
set(gca,'fontsize',14)
title('Number of Trials To Reach Threshold','fontsize',20)
xlabel('N','fontsize',14)
ylabel('Number of Trials','fontsize',14)
ylim([0 250])

end

function trls = calcStairThresh(s0,ds,ntrls,M,N,error)

s = zeros(ntrls,1);
s(1) = s0;

for i = 1:ntrls
    p = erf(s(i)/sqrt(2))/2 + 0.5;
    s(i+1) = s(i) + ds*(1-p)^M - ds*p^N;
end

rootvec = zeros(1,N+1);
rootvec(end-N) = 1;
for i = 0:M
    rootvec(end-i) = rootvec(end-i) + ...
        (-1)^(i+1)*factorial(M)/(factorial(i)*factorial(M-i));
end
allroots = roots(rootvec);
accuracy = real(allroots(abs(imag(allroots)) < 0.001 & allroots >= 0.499 & allroots < 1));

curve = erf(s./sqrt(2))/2 + 0.5;
threshpct = (1 + error)*accuracy;
%threshpct = accuracy + error;

trls = find(curve <= threshpct); trls = trls(1)-1;

end

function trls = BlockStairSim(s0,ds,ntrls,M,N,error)

s = zeros(ntrls,1);
s(1) = s0;

for i = 1:ntrls
    p = erf(s(i)/sqrt(2))/2 + 0.5; % added missing 0.5, possible bug int he previous version??
    s(i+1) = s(i) + ds*(p*(1-p)^M/(1-(1-p)^M)) - ds*((1-p)*p^N/(1-p^N));
end

rootvec = zeros(1,N+M);
rootvec(end-(N-1)) = 1;
for i = 0:M
    rootvec(end-(i+N-1)) = rootvec(end-(i+N-1)) + ...
        (-1)^(i+1)*factorial(M)/(factorial(i)*factorial(M-i));
end
for i = 0:M-1
    rootvec(end-i) = rootvec(end-i) + ...
        (-1)^(i+1)*factorial(M-1)/(factorial(i)*factorial(M-1-i));
    rootvec(end-(i+N)) = rootvec(end-(i+N)) + ...
        (-1)^(i)*factorial(M-1)/(factorial(i)*factorial(M-1-i));
end
allroots = roots(rootvec);
accuracy = real(allroots(abs(imag(allroots)) < 0.001 & allroots >= 0.499 & allroots < 1));

curve = erf(s./sqrt(2))/2 + 0.5; % added missing 0.5, possible bug int he previous version??
threshpct = (1 + error)*accuracy;
%threshpct = accuracy + error;

trls = find(curve <= threshpct); trls = trls(1)-1;

end

function trls = NonConsecStairSim(s0,ds,ntrls,M,N,error)

s = zeros(ntrls,1);
s(1) = s0;

for i = 1:ntrls
    p = erf(s(i)/sqrt(2))/2 + 0.5; % added missing 0.5, possible bug int he previous version??
    s(i+1) = s(i) + (1-p)*ds/M - p*ds/N;
end

accuracy = N/(M+N);

curve = erf(s./sqrt(2))/2 + 0.5; % added missing 0.5, possible bug int he previous version??
threshpct = (1 + error)*accuracy;
%threshpct = accuacy + error;

trls = find(curve <= threshpct); trls = trls(1)-1;

end

function trls = WeightedStairSim(s0,ds,ntrls,gamma,error)

s = zeros(ntrls,1);
s(1) = s0;
% ds = ds/((gamma + 1)/2); % preserve arithmetic mean
% ds = ds/sqrt(gamma*1); % preserve geo mean

for i = 1:ntrls
    p = erf(s(i)/sqrt(2))/2 + 0.5; % added missing 0.5, possible bug int he previous version??
    s(i+1) = s(i) + (1-p)*ds/1 - p*ds/gamma;
end

accuracy = gamma/(gamma+1);

curve = erf(s./sqrt(2))/2 + 0.5; % added missing 0.5, possible bug int he previous version??
threshpct = (1 + error)*accuracy;
%threshpct = accuracy + error;
if threshpct<0.5, threshpct = 0.5; end

trls = find(curve <= threshpct); trls = trls(1)-1;

end