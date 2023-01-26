close all
clear

set(groot,'defaultLineLineWidth',2);
set(groot,'defaultFigureColor','w');
set(groot,'defaultTextFontsize',22);
set(groot,'defaultAxesFontsize',22);
set(groot,'defaultPolarAxesFontsize',22);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultPolarAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesLineWidth',1);
%% load data
load('evasionData.mat');
phi = data.Phi(isfinite(data.U));
u   = data.U(isfinite(data.U));
v   = data.V(isfinite(data.U))*0.01;
lam = data.Lambda(isfinite(data.U));
theta = data.Delta(isfinite(data.U));

theta = wrapTo2Pi(theta);
%% S2

load('s2_opt_result.mat');
theta_s2_d = phi + sign(pi-phi)*pi;

theta_s2_d= wrapTo2Pi(theta_s2_d);
sig_phi = etaopt(1);
sig_theta = etaopt(2);
% sig_phi = 0.001;
% sig_theta = 0.001;

r = 500;
N = length(theta);
theta_s2 = zeros(N,r);
for j = 1:r
for i = 1:N
    %     phih = phi(i);
    phih = randraw('vonmises', [phi(i),1/sig_phi^2], [1 1]);
    thetah = phih + sign(pi-phih)*pi;
    theta_s2(i,j) = wrapTo2Pi(randraw('vonmises', [thetah,1/sig_theta^2], [1 1]));
    %     theta_s2(i) = wrapTo2Pi(thetah);
end
end
theta_s2 = wrapTo2Pi(theta_s2);
% figure, polarhistogram(theta_s2(:,1),36,'Normalization','pdf');

%% S3

load('s3_opt_result.mat')
sig_phi = etaopt(1);
sig_lam = etaopt(2);
sig_theta = etaopt(3);

theta_s3 = zeros(N,r);
for j = 1:r
parfor i = 1:N
    phih = randraw('vonmises', [phi(i),1/sig_phi^2], [1,1]);
    lamh = randraw('vonmises', [lam(i),1/sig_lam^2], [1,1]);
    %     phih = phi(i);
    %     lamh = lam(i);
    thetah = phih + lamh + pi - sign(lamh)*pi/2;
    theta_s3(i,j) = wrapTo2Pi(randraw('vonmises', [thetah,1/sig_theta^2], [1,1]));
    %     theta_s3(i,r) = wrapTo2Pi(thetah);
end
end

% figure, polarhistogram(theta_s3,36,'Normalization','pdf');
%%
k = 1;
X = [theta;phi:lam];
for l = 1:r
    Y = [theta_s2(:,l)'; phi; lam];
Y2 = theta_s2(:,l)';
% Y3 = theta_s3(:,l)';
kl_s2(l) = KLDivergence(X,Y2,k);
% kl_s3(l) = KLDivergence(X,Y3,k);
end
Y_d = [theta_s2_d;phi;lam];
kl_s2(l) = KLDivergence(X,Y2,k);
kl_s2_d = KLDivergence(X,Y_d,k);
%%
[h2,p2,ci2,stats2] = ttest(kl_s2);
[h3,p3,ci3,stats3] = ttest(kl_s3);
[h23,p23,ci23,stats23] = ttest(kl_s2,kl_s3);
