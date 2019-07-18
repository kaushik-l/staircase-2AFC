function NumRealRoots = ComputeNumRealRoots(polycoeffs,interval)

a = interval(1); 
b = interval(2);

%% compute sequence of polynomials P_i
P{1} = polycoeffs(:); N(1) = numel(P{1})-1;
P{2} = flip(0:N(1))'.*P{1}(:); P{2} = P{2}(1:end-1); N(2) = numel(P{2})-1;
while N(end)>0
    [~,R] = deconv(P{end-1},P{end}); P{end+1} = -R(find(R~=0,1):end);
    N(end+1) = numel(P{end})-1;
end

%% count number of sign variations
for i=1:numel(P)
    P_a(i) = sum(flip(a.^(0:N(i)))'.*P{i});
    P_b(i) = sum(flip(b.^(0:N(i)))'.*P{i});
end
V_a = floor(sum(abs(diff(sign(P_a))))/2);
V_b = floor(sum(abs(diff(sign(P_b))))/2);
% V_a = sum(abs(diff(sign(P_a)))==2);
% V_b = sum(abs(diff(sign(P_b)))==2);
% V_a = sum(abs(diff(sign(P_a)))>0);
% V_b = sum(abs(diff(sign(P_b)))>0);
NumRealRoots = V_a - V_b;