clear;
close all;

set(groot,'defaultLineLineWidth',2);
set(groot,'defaultFigureColor','w');
set(groot,'defaultTextFontsize',18);
set(groot,'defaultAxesFontsize',18);
set(groot,'defaultPolarAxesFontsize',8);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultPolarAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesLineWidth',1);

load('evasionData.mat');
phi = data.Phi(isfinite(data.U));
u   = data.U(isfinite(data.U));
v   = data.V(isfinite(data.U))*0.01;
lam = data.Lambda(isfinite(data.U));
theta = data.Delta(isfinite(data.U));

theta = wrapToPi(theta);
%% Contralateral
theta_CL = sign(phi-pi)*pi/2;

theta_CL= wrapToPi(theta_CL);

%% Antipodal
theta_AP = phi + sign(pi-phi)*pi;

theta_AP= wrapToPi(theta_AP);
%% Parallel
theta_PL = phi + lam + pi;

theta_PL= wrapToPi(theta_PL);
%% Distance-optimal, u/v = 0.5
chi = acos(0.5);
theta_DO = phi + lam + pi - ((lam>=0)-0.5)*chi*2;

theta_DO= wrapToPi(theta_DO);
%% Orthogonal
chi = acos(0);
theta_OT = phi + lam + pi - ((lam>=0)-0.5)*chi*2;

theta_OT= wrapToPi(theta_OT);
%% Distance-optimal, u/v from data (real)
% chi = acos(0.1);
% theta_s4_3 = phi + lam + pi - ((lam>=0)-0.5)*chi*2;
chi = real(acos(u./v));
theta_DOR = phi + lam + pi - sign(lam).*chi;
theta_DOR= wrapTo2Pi(theta_DOR);
%% saperation of data according to different v
phi_02 = phi(v == 0.02);
phi_11 = phi(v == 0.11);
phi_20 = phi(v == 0.2);
lam_02 = lam(v == 0.02);
lam_11 = lam(v == 0.11);
lam_20 = lam(v == 0.2);
u_02 = u(v == 0.02);
u_11 = u(v == 0.11);
u_20 = u(v == 0.2);
theta_02 = theta(v == 0.02);
theta_11 = theta(v == 0.11);
theta_20 = theta(v == 0.2);

theta_CL_02 = theta_CL(v == 0.02);
theta_CL_11 = theta_CL(v == 0.11);
theta_CL_20 = theta_CL(v == 0.2);
theta_AP_02 = theta_AP(v == 0.02);
theta_AP_11 = theta_AP(v == 0.11);
theta_AP_20 = theta_AP(v == 0.2);
theta_PL_02 = theta_PL(v == 0.02);
theta_PL_11 = theta_PL(v == 0.11);
theta_PL_20 = theta_PL(v == 0.2);
theta_DO_02 = theta_DO(v == 0.02);
theta_DO_11 = theta_DO(v == 0.11);
theta_DO_20 = theta_DO(v == 0.2);
theta_OT_02 = theta_OT(v == 0.02);
theta_OT_11 = theta_OT(v == 0.11);
theta_OT_20 = theta_OT(v == 0.2);
theta_DOR_02 = theta_DOR(v == 0.02);
theta_DOR_11 = theta_DOR(v == 0.11);
theta_DOR_20 = theta_DOR(v == 0.2);
%% all
fig = figure(1);
fig.Color = 'w';
fig.Position = [166.6000 281.8000 500 115];
ax_CL = subplot(1,6,1,polaraxes);
ph_CL = polarhistogram(theta_CL,36);


ax_AP = subplot(1,6,2,polaraxes);
ph_AP = polarhistogram(theta_AP,36);

ax_PL = subplot(1,6,3,polaraxes);
ph_PL = polarhistogram(theta_PL,36);

ax_DO = subplot(1,6,4,polaraxes);
ph_DO = polarhistogram(theta_DO,36);

ax_OT = subplot(1,6,5,polaraxes);
ph_OT = polarhistogram(theta_OT,36);

ax_exp = subplot(1,6,6,polaraxes);
ph_exp = polarhistogram(theta,36);

ax_CL.FontSize = 8;
ax_AP.FontSize = 8;
ax_PL.FontSize = 8;
ax_DO.FontSize = 8;
ax_OT.FontSize = 8;
title(ax_CL,'Contralateral');
title(ax_AP,'Antipodal');
title(ax_PL,'Parallel');
title(ax_DO,'Distance-optimal (u/v =  0.5)');
title(ax_OT,'Orthogonal');
title(ax_exp,'Experiment');

rticklabels(ax_CL,{});
rticklabels(ax_AP,{});
rticklabels(ax_PL,{});
rticklabels(ax_DO,{});
rticklabels(ax_OT,{});
rticklabels(ax_exp,{});
ax_CL.RGrid = 'off';
ax_AP.RGrid = 'off';
ax_PL.RGrid = 'off';
ax_DO.RGrid = 'off';
ax_OT.RGrid = 'off';
ax_exp.RGrid = 'off';
ax_CL.ThetaTick = [0 45 90 135 180 225 270 315];
ax_CL.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
ax_AP.ThetaTick = [0 45 90 135 180 225 270 315];
ax_AP.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
ax_PL.ThetaTick = [0 45 90 135 180 225 270 315];
ax_PL.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
ax_DO.ThetaTick = [0 45 90 135 180 225 270 315];
ax_DO.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
ax_OT.ThetaTick = [0 45 90 135 180 225 270 315];
ax_OT.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
ax_exp.ThetaTick = [0 45 90 135 180 225 270 315];
ax_exp.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
ph_CL.FaceColor = [0.2 0.2 0.2];
ph_AP.FaceColor = [0.2 0.2 0.2];
ph_PL.FaceColor = [0.2 0.2 0.2];
ph_DO.FaceColor = [0.2 0.2 0.2];
ph_OT.FaceColor = [0.2 0.2 0.2];
ph_exp.FaceColor = [0.2 0.2 0.2];

%% 3-D KL-Divergence per v
k = 5;
X = [theta_02; phi_02; lam_02]; % s1
Y = [theta_CL_02; phi_02; lam_02]; % s1
kl_CL(1) = KLDivergence_angle(X,Y,k);
Y = [theta_AP_02; phi_02; lam_02]; % s2
kl_AP(1) = KLDivergence_angle(X,Y,k);
Y = [theta_PL_02; phi_02; lam_02]; % s3
kl_PL(1) = KLDivergence_angle(X,Y,k);
Y = [theta_DO_02; phi_02; lam_02]; % s4-1
kl_DO(1) = KLDivergence_angle(X,Y,k);
Y = [theta_OT_02; phi_02; lam_02]; % s4-2
kl_OT(1) = KLDivergence_angle(X,Y,k);
% Y = [theta_DOR_02; phi_02; lam_02]; % s4-3
% kl_DOR(1) = KLDivergence(X,Y,k);

X = [theta_11; phi_11; lam_11]; % s1
Y = [theta_CL_11; phi_11; lam_11]; % s1
kl_CL(2) = KLDivergence_angle(X,Y,k);
Y = [theta_AP_11; phi_11; lam_11]; % s2
kl_AP(2) = KLDivergence_angle(X,Y,k);
Y = [theta_PL_11; phi_11; lam_11]; % s3
kl_PL(2) = KLDivergence_angle(X,Y,k);
Y = [theta_DO_11; phi_11; lam_11]; % s4-1
kl_DO(2) = KLDivergence_angle(X,Y,k);
Y = [theta_OT_11; phi_11; lam_11]; % s4-2
kl_OT(2) = KLDivergence_angle(X,Y,k);
% Y = [theta_DOR_11; phi_11; lam_11]; % s4-3
% kl_DOR(2) = KLDivergence(X,Y,k);

X = [theta_20; phi_20; lam_20]; % s1
Y = [theta_CL_20; phi_20; lam_20]; % s1
kl_CL(3) = KLDivergence_angle(X,Y,k);
Y = [theta_AP_20; phi_20; lam_20]; % s2
kl_AP(3) = KLDivergence_angle(X,Y,k);
Y = [theta_PL_20; phi_20; lam_20]; % s3
kl_PL(3) = KLDivergence_angle(X,Y,k);
Y = [theta_DO_20; phi_20; lam_20]; % s4-1
kl_DO(3) = KLDivergence_angle(X,Y,k);
Y = [theta_OT_20; phi_20; lam_20]; % s4-2
kl_OT(3) = KLDivergence_angle(X,Y,k);
% Y = [theta_DOR_20; phi_20; lam_20]; % s4-3
% kl_DOR(3) = KLDivergence(X,Y,k);



X = [theta; phi; lam]; % s1
Y = [theta_CL; phi; lam]; % s1
kl_CL(4) = KLDivergence_angle(X,Y,k);
Y = [theta_AP; phi; lam]; % s2
kl_AP(4) = KLDivergence_angle(X,Y,k);
Y = [theta_PL; phi; lam]; % s3
kl_PL(4) = KLDivergence_angle(X,Y,k);
Y = [theta_DO; phi; lam]; % s4-1
kl_DO(4) = KLDivergence_angle(X,Y,k);
Y = [theta_OT; phi; lam]; % s4-2
kl_OT(4) = KLDivergence_angle(X,Y,k);
% Y = [theta_DOR; phi; lam]; % s4-3
% kl_DOR(4) = KLDivergence(X,Y,k);

K = [kl_CL;kl_AP;kl_PL;kl_DO;kl_OT]';
fKL = figure;bp = bar([K zeros(4,2)]); fKL.Position = [325 279 390 156];box off;
bp(1).FaceColor = [0.25,0.6,0.9];bp(1).EdgeColor = 'none';
bp(2).FaceColor = [1,0.7,0.3];bp(2).EdgeColor = 'none';
bp(3).FaceColor = [0.1,0.7,0.25];bp(3).EdgeColor = 'none';
bp(4).FaceColor = 'k';bp(4).EdgeColor = 'none';
% legend('dataset 1','dataset 2','dataset 3','all data');
% xticks('');
% title(['k = ' num2str(k)]);
%% finding best chi
k = 20;
figure();hold on;
chi = (0:18)*pi/36;
kl = zeros(1,19);
kl_02 = zeros(1,19);
kl_11 = zeros(1,19);
kl_20 = zeros(1,19);
for j = 1:19
    theta_s4 = phi + lam + pi - ((lam>=0)-0.5)*chi(j)*2;
    theta_s4_02 = phi_02 + lam_02 + pi - ((lam_02>=0)-0.5)*chi(j)*2;
    %     theta_s4_02= wrapTo2Pi(theta_s4_02);
    theta_s4_11 = phi_11 + lam_11 + pi - ((lam_11>=0)-0.5)*chi(j)*2;
    %     theta_s4_11= wrapTo2Pi(theta_s4_11);
    theta_s4_20 = phi_20 + lam_20 + pi - ((lam_20>=0)-0.5)*chi(j)*2;
    %     theta_s4_20= wrapTo2Pi(theta_s4_20);
    X = [theta; phi; lam];
    Y = [theta_s4; phi; lam];
    X_02 = [theta_02; phi_02; lam_02];
    X_11 = [theta_11; phi_11; lam_11];
    X_20 = [theta_20; phi_20; lam_20];
    Y_02 = [theta_s4_02; phi_02; lam_02];
    Y_11 = [theta_s4_11; phi_11; lam_11];
    Y_20 = [theta_s4_20; phi_20; lam_20];
    kl(j) = KLDivergence(X,Y,k);
    kl_02(j) = KLDivergence(X_02,Y_02,k);
    kl_11(j) = KLDivergence(X_11,Y_11,k);
    kl_20(j) = KLDivergence(X_20,Y_20,k);
end
plot(chi,kl,'ko-');
plot(chi,kl_02,'ro-');
plot(chi,kl_11,'go-');
plot(chi,kl_20,'bo-');
xlabel('$\chi$ [rad]');
xticks([0,pi/4,pi/2]); xticklabels({'0','$\pi/4$','$\pi/2$'});
ylabel('KL divergence');
title(['k = ' num2str(k)]);
legend('all','slow predator','medium speed predator','fast predator');
%% 4-D KL-Divergence Bootstrapping
% k = 1;
% 
% kl_his = zeros(15,4);
% for j = 1:200
%     r = randperm(699,600);
%     X = [theta(r); phi(r); lam(r); v(r)]; % s1
%     Y = [theta_s1(r); phi(r); lam(r); v(r)]; % s1
%     % X = theta(r); % s1
%     % Y = theta_s1(r); % s1
%     
%     kl_s1 = KLDivergence(X,Y,k);
%     kl_his(j,1) = kl_s1;
%     
%     r = randperm(699,600);
%     X = [theta(r); phi(r); lam(r); v(r)];
%     Y = [theta_s2(r); phi(r); lam(r); v(r)]; % s2
%     % X = theta(r); % s2
%     % Y = theta_s2(r); % s2
%     kl_s2 = KLDivergence(X,Y,k);
%     kl_his(j,2) = kl_s2;
%     
%     r = randperm(699,600);
%     X = [theta(r); phi(r); lam(r); v(r)];
%     Y = [theta_s3(r); phi(r); lam(r); v(r)]; % s3
%     % X = theta(r); % s3
%     % Y = theta_s3(r); % s3
%     kl_s3 = KLDivergence(X,Y,k);
%     kl_his(j,3) = kl_s3;
%     
%     r = randperm(699,600);
%     X = [theta(r); phi(r); lam(r); v(r)];
%     Y = [theta_s4(r); phi(r); lam(r); v(r)]; % s4
%     % X = theta(r); % s4
%     % Y = theta_s4(r); % s4
%     
%     kl_s4 = KLDivergence(X,Y,k);
%     kl_his(j,4) = kl_s4;
% end



%%
figure,histogram(kl_his(:,1),10);
hold on;
histogram(kl_his(:,2),10);
histogram(kl_his(:,3),10);
histogram(kl_his(:,4),10);
legend('1','2','3','4');
%%
% theta = theta(randperm(699));
% for j = 1:200
%     X = theta(1:350);
%     Y = theta(351:end);
% r = randperm(349,300);
% X = X(r);
% % r = randperm(350,300);
% Y = Y(r);
% kl_test(j) = KLDivergence(X,Y,k);
% end
%%
[h1,p1,ci1,stats1] = ttest(kl_his(:,1));
[h2,p2,ci2,stats2] = ttest(kl_his(:,2));
[h3,p3,ci3,stats3] = ttest(kl_his(:,3));
[h4,p4,ci4,stats4] = ttest(kl_his(:,4));
%%
[h23,p23,ci23,stats23] = ttest(kl_his(:,2),kl_his(:,3));
[h24,p24,ci24,stats24] = ttest(kl_his(:,2),kl_his(:,4));
[h34,p34,ci34,stats34] = ttest(kl_his(:,3),kl_his(:,4));

