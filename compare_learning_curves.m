function compare_learning_curves(s0,ds,ntrls,M,N)

[rel_curve_cont,rel_int_curve_cont] = calcStairThresh(s0,ds,ntrls,M,N);
[rel_curve_block,rel_int_curve_block] = BlockStairSim(s0,ds,ntrls,M,N);
[rel_curve_nc,rel_int_curve_nc] = NonConsecStairSim(s0,ds,ntrls,M,N);

figure, hold on
plot(rel_curve_cont), plot(rel_curve_block), plot(rel_curve_nc)
set(gca,'fontsize',14)
title(['M = ',num2str(M),', N = ',num2str(N)])
xlabel('Trials','fontsize',14)
ylabel('Relative Stimulus Magnitude','fontsize',14)
legend('Continuous','Blocked','Non-Consecutive','location','southeast')
ylim([0 1])

% figure, hold on
% plot(rel_int_curve_cont), plot(rel_int_curve_block), plot(rel_int_curve_nc)
% set(gca,'fontsize',14)
% title(['Integral, M = ',num2str(M),', N = ',num2str(N)])
% xlabel('Trials','fontsize',14)
% ylabel('Cumulative Area/Total Area','fontsize',14)
% legend('Continuous','Blocked','Non-Consecutive','location','southeast')
% ylim([0 1])

end

function [rel_curve,rel_int_curve] = calcStairThresh(s0,ds,ntrls,M,N)

s = zeros(ntrls,1);
s(1) = s0;

for i = 1:ntrls
    p = erf(s(i)/sqrt(2));
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
threshVal = sqrt(2)*erfinv(accuracy);

integral = trapz(s-threshVal);

rel_curve = 1 - (s - threshVal)./(s0 - threshVal);
rel_int_curve = cumtrapz(s-threshVal)./integral;

end

function [rel_curve,rel_int_curve] = BlockStairSim(s0,ds,ntrls,M,N)

s = zeros(ntrls,1);
s(1) = s0;

for i = 1:ntrls
    p = erf(s(i)/sqrt(2));
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
threshVal = sqrt(2)*erfinv(accuracy);

integral = trapz(s-threshVal);

rel_curve = 1 - (s - threshVal)./(s0 - threshVal);
rel_int_curve = cumtrapz(s-threshVal)./integral;

end

function [rel_curve,rel_int_curve] = NonConsecStairSim(s0,ds,ntrls,M,N)

s = zeros(ntrls,1);
s(1) = s0;

for i = 1:ntrls
    p = erf(s(i)/sqrt(2));
    s(i+1) = s(i) + (1-p)*ds/M - p*ds/N;
end

accuracy = N/(M+N);
threshVal = sqrt(2)*erfinv(accuracy);

integral = trapz(s-threshVal);

rel_curve = 1 - (s - threshVal)./(s0 - threshVal);
rel_int_curve = cumtrapz(s-threshVal)./integral;

end