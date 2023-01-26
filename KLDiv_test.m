close all;
clear;
%% test of KL divergence estimation
N = 1000;
Kbound = 50;
v1 = normrnd(0 , 1 , 1 , N);
v2 = normrnd(0 , 1 , 1 , N);
v3 = normrnd(0 , 1.1 , 1 , N);
v4 = normrnd(0.1 , 1 , 1 , N);
kl1 = zeros(1,Kbound);
kl2 = zeros(1,Kbound);
kl3 = zeros(1,Kbound);
for k = 1:Kbound
kl1(k) = KLDivergence(v1,v2,k);
kl2(k) = KLDivergence(v1,v3,k);
kl3(k) = KLDivergence(v1,v4,k);
end
figure,hold on;
plot(1:Kbound,kl1,'ro-');
plot(1:Kbound,kl2,'bo-');
plot(1:Kbound,kl3,'go-');
kl1_sol = klNorm(0 , 1 , 0 , 1);
kl2_sol = klNorm(0 , 1 , 0 , 1.1);
kl3_sol = klNorm(0 , 1 , 0.1 , 1);
yline(kl1_sol,'r--','linewidth',1.5);
yline(kl2_sol,'b--','linewidth',1.5);
yline(kl3_sol,'g--','linewidth',1.5);
legend('same','diff var','diff mean');
xlabel('choice of k');
ylabel('KL divergence estimate');
kl4_sol = klNorm(0 , 1 , 0 , 2);
%%
function sol = klNorm(mu1,sig1,mu2,sig2)
sol = log(sig2/sig1) + 1/2*((sig1/sig2)^2 - 1 + (mu1-mu2)^2/sig2^2);
end
% sol = log(sig2/sig1) + (sig1^2+(mu1-mu2)^2)/sig2^2/2 - 1/2;