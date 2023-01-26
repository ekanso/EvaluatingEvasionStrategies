clear;
close all;
clc;

set(groot,'defaultLineLineWidth',2);
set(groot,'defaultFigureColor','w');
set(groot,'defaultTextFontsize',10);
set(groot,'defaultAxesFontsize',10);
set(groot,'defaultPolarAxesFontsize',10);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultPolarAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesLineWidth',1);

load('evasionData.mat');
phi = wrapTo2Pi(data.Phi(isfinite(data.U)));
u   = data.U(isfinite(data.U));
v   = data.V(isfinite(data.U))*0.01;
lam = data.Lambda(isfinite(data.U));
theta = data.Delta(isfinite(data.U));


theta = wrapToPi(theta);
%% S1
theta_CL = sign(phi-pi)*pi/2;

theta_CL = wrapToPi(theta_CL);
% theta_s1(theta_s1>2*pi) = theta_s1(theta_s1>2*pi) - 2*pi;
% theta_s1(theta_s1<0) = theta_s1(theta_s1<0) + 2*pi;

%% S2
theta_AP = phi + sign(pi-phi)*pi;

% theta_s2(theta_s2>2*pi) = theta_s2(theta_s2>2*pi) - 2*pi;
% theta_s2(theta_s2<0) = theta_s2(theta_s2<0) + 2*pi;
theta_AP = wrapToPi(theta_AP);
%% S3
theta_OT = phi + lam + pi - sign(lam)*pi/2;

% theta_s3(theta_s3>2*pi) = theta_s3(theta_s3>2*pi) - 2*pi;
% theta_s3(theta_s3<0) = theta_s3(theta_s3<0) + 2*pi;
theta_OT = wrapToPi(theta_OT);
%% S4
% u_mod = u;
% u_mod(u>v) = v(u>v);
chi = real(acos(u./v));

theta_DOR = phi + lam + pi - sign(lam).*chi;

% theta_s4(theta_s4>2*pi) = theta_s4(theta_s4>2*pi) - 2*pi;
% theta_s4(theta_s4<0) = theta_s4(theta_s4<0) + 2*pi;
theta_DOR = wrapToPi(theta_DOR);
%% Distance-optimal
% u_mod = u;
% u_mod(u>v) = v(u>v);
chi = acos(0.5);

theta_DO = phi + lam + pi - sign(lam).*chi;

% theta_s4(theta_s4>2*pi) = theta_s4(theta_s4>2*pi) - 2*pi;
% theta_s4(theta_s4<0) = theta_s4(theta_s4<0) + 2*pi;
theta_DO = wrapToPi(theta_DO);
%% saperation of data according to different speed ratio
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
theta_OT_02 = theta_OT(v == 0.02);
theta_OT_11 = theta_OT(v == 0.11);
theta_OT_20 = theta_OT(v == 0.2);
theta_DOR_02 = theta_DOR(v == 0.02);
theta_DOR_11 = theta_DOR(v == 0.11);
theta_DOR_20 = theta_DOR(v == 0.2);

%% correlation test between phi and lambda
figure('position',[440 496 450 300]);
subplot(1,3,1),plot(abs(phi_02),abs(lam_02),'.');
axis equal;xlim([0,pi]);ylim([0,pi]);
xticks([0,pi/2,pi]);yticks([0,pi/2,pi]);
xticklabels({'0','$\pi$/2','$\pi$'});yticklabels({'0','$\pi$/2','$\pi$'});
title('slow predator','fontSize',20); xlabel('abs$(\phi)$','fontSize',10); ylabel('abs$(\lambda)$','fontSize',10);
subplot(1,3,2),plot(abs(phi_11),abs(lam_11),'.');
axis equal;xlim([0,pi]);ylim([0,pi]);xticks([0,pi/2,pi]);yticks([0,pi/2,pi]);
xticklabels({'0','$\pi$/2','$\pi$'});yticklabels({'0','$\pi$/2','$\pi$'});
title('mid-speed predator','fontSize',20); xlabel('abs$(\phi)$','fontSize',10); ylabel('abs$(\lambda)$','fontSize',10);
subplot(1,3,3),plot(abs(phi_20),abs(lam_20),'.');
axis equal;xlim([0,pi]);ylim([0,pi]);xticks([0,pi/2,pi]);yticks([0,pi/2,pi]);
xticklabels({'0','$\pi$/2','$\pi$'});yticklabels({'0','$\pi$/2','$\pi$'});
title('fast predator','fontSize',20); xlabel('abs$(\phi)$','fontSize',10); ylabel('abs$(\lambda)$','fontSize',10);

PL_total = [abs(phi)',abs(lam)'];
PL_slow = [abs(phi_02)',abs(lam_02)'];
PL_mid = [abs(phi_11)',abs(lam_11)'];
PL_fast = [abs(phi_20)',abs(lam_20)'];

PT_total = [wrapTo2Pi(phi)',theta'];
PT_slow = [wrapTo2Pi(phi_02)',theta_02'];
PT_mid = [wrapTo2Pi(phi_11)',theta_11'];
PT_fast = [wrapTo2Pi(phi_20)',theta_20'];
%%

% [rhoPL pvalPL] = circ_corrcc(abs(phi)',abs(lam)');
% [rhoPLs pvalPLs] = circ_corrcc(abs(phi_02)',abs(lam_02)');
% [rhoPLm pvalPLm] = circ_corrcc(abs(phi_11)',abs(lam_11)');
% [rhoPLf pvalPLf] = circ_corrcc(abs(phi_20)',abs(lam_20)');

[rhoPL pvalPL] = circ_corrcc((phi)',(lam)');
[rhoPLs pvalPLs] = circ_corrcc((phi_02)',(lam_02)');
[rhoPLm pvalPLm] = circ_corrcc((phi_11)',(lam_11)');
[rhoPLf pvalPLf] = circ_corrcc((phi_20)',(lam_20)');

[rhoPT pvalPT] = circ_corrcc(phi',theta');
[rhoPTs pvalPTs] = circ_corrcc(phi_02',theta_02');
[rhoPTm pvalPTm] = circ_corrcc(phi_11',theta_11');
[rhoPTf pvalPTf] = circ_corrcc(phi_20',theta_20');

[rhoLT pvalLT] = circ_corrcc(lam',theta');
[rhoLTs pvalLTs] = circ_corrcc(lam_02',theta_02');
[rhoLTm pvalLTm] = circ_corrcc(lam_11',theta_11');
[rhoLTf pvalLTf] = circ_corrcc(lam_20',theta_20');


psi = phi+lam;
psi_02 = phi_02 + lam_02;
psi_11 = phi_11 + lam_11;
psi_20 = phi_20 + lam_20;

[rhoPLT, pvalPLT] = circ_corrcc(psi',theta');
[rhoPLTs, pvalPLTs] = circ_corrcc(psi_02',theta_02');
[rhoPLTm, pvalPLTm] = circ_corrcc(psi_11',theta_11');
[rhoPLTf, pvalPLTf] = circ_corrcc(psi_20',theta_20');

[rhoPLT_pL, pvalPLT_pL] = circ_corrcc(psi(lam>0)',theta(lam>0)');
[rhoPLTs_pL, pvalPLTs_pL] = circ_corrcc(psi_02(lam_02>0)',theta_02(lam_02>0)');
[rhoPLTm_pL, pvalPLTm_pL] = circ_corrcc(psi_11(lam_11>0)',theta_11(lam_11>0)');
[rhoPLTf_pL, pvalPLTf_pL] = circ_corrcc(psi_20(lam_20>0)',theta_20(lam_20>0)');

[rhoPLT_nL, pvalPT_nL] = circ_corrcc(psi(lam<0)',theta(lam<0)');
[rhoPLTs_nL, pvalPLTs_nL] = circ_corrcc(psi_02(lam_02<0)',theta_02(lam_02<0)');
[rhoPLTm_nL, pvalPLTm_nL] = circ_corrcc(psi_11(lam_11<0)',theta_11(lam_11<0)');
[rhoPLTf_nL, pvalPLTf_nL] = circ_corrcc(psi_20(lam_20<0)',theta_20(lam_20<0)');
disp('phi-lam');
[rhoPL,rhoPLs,rhoPLm,rhoPLf]
[pvalPL,pvalPLs,pvalPLm,pvalPLf]
disp('phi-theta');
[rhoPT,rhoPTs,rhoPTm,rhoPTf]
[pvalPT,pvalPTs,pvalPTm,pvalPTf]
disp('lam-theta');
[rhoLT,rhoLTs,rhoLTm,rhoLTf]
[pvalLT,pvalLTs,pvalLTm,pvalLTf]
disp('phi+lam,psi-theta');
[rhoPLT,rhoPLTs,rhoPLTm,rhoPLTf]
[pvalPLT,pvalPLTs,pvalPLTm,pvalPLTf]
disp('phi+lam,psi-theta(lam>0)');
[rhoPLT_pL,rhoPLTs_pL,rhoPLTm_pL,rhoPLTf_pL]
[pvalPLT_pL,pvalPLTs_pL,pvalPLTm_pL,pvalPLTf_pL]
disp('phi+lam,psi-theta(lam<0)');
[rhoPLT_nL,rhoPLTs_nL,rhoPLTm_nL,rhoPLTf_nL]
[pvalPT_nL,pvalPLTs_nL,pvalPLTm_nL,pvalPLTf_nL]

% PL_total = [abs(phi)',abs(lam)'];
% PL_slow = [abs(phi_02)',abs(lam_02)'];
% PL_mid = [abs(phi_11)',abs(lam_11)'];
% PL_fast = [abs(phi_20)',abs(lam_20)'];
% B_slow = [abs(wrapTo2Pi(phi_02+lam_02))',abs(theta_02)'];
% B_mid = [abs(wrapTo2Pi(phi_11+lam_11))',abs(theta_11)'];
% B_fast = [abs(wrapTo2Pi(phi_20+lam_20))',abs(theta_20)'];
% 
% B_slowP = [wrapTo2Pi(phi_02(lam_02>0)+lam_02(lam_02>0))',theta_02(lam_02>0)'];
% B_midP = [wrapTo2Pi(phi_11(lam_11>0)+lam_11(lam_11>0))',theta_11(lam_11>0)'];
% B_fastP = [wrapTo2Pi(phi_20(lam_20>0)+lam_20(lam_20>0))',theta_20(lam_20>0)'];
% B_slowM = [wrapTo2Pi(phi_02(lam_02<0)+lam_02(lam_02<0))',theta_02(lam_02<0)'];
% B_midM = [wrapTo2Pi(phi_11(lam_11<0)+lam_11(lam_11<0))',theta_11(lam_11<0)'];
% B_fastM = [wrapTo2Pi(phi_20(lam_20<0)+lam_20(lam_20<0))',theta_20(lam_20<0)'];

% [R_t,P_t] = corrcoef(PL_total)
% [R_s,P_s] = corrcoef(PL_slow)
% % RS_s = corr(A_slow,'Type','Spearman')
% [R_m,P_m] = corrcoef(PL_mid)
% % RS_m = corr(A_mid,'Type','Spearman')
% [R_f,P_f] = corrcoef(PL_fast)
% % RS_f = corr(A_fast,'Type','Spearman')

% [RPT_t,PPT_t] = corrcoef(PT_total)
% [RPT_s,PPT_s] = corrcoef(PT_slow)
% [RPT_m,PPT_m] = corrcoef(PT_mid)
% [RPT_f,PPT_f] = corrcoef(PT_fast)
% 
% [RP_s,PP_s] = corrcoef(B_slowP)
% [RP_m,PP_m] = corrcoef(B_midP)
% [RP_f,PP_f] = corrcoef(B_fastP)
% 
% [RM_s,PM_s] = corrcoef(B_slowM)
% [RM_m,PM_m] = corrcoef(B_midM)
% [RM_f,PM_f] = corrcoef(B_fastM)
%%

%% AP vs Experiment
figure,plot(phi_20,theta_20,'b.','markersize',3);
hold on;
plot(phi_20,theta_AP_20,'rd','markersize',3);
axis equal;
xlim([0,2*pi]);
ylim([-pi,pi]);
xticks([0,pi,2*pi]);
yticks([-pi,0,pi]);

%% polar histograms for all datasets
% fig = figure(1);
% fig.Color = 'w';
% fig.Position = [166.6000 281.8000 1.7008e+03 503.2000];
% ax_s1 = subplot(1,5,1,polaraxes);
% ph_1 = polarhistogram(theta_s1,36);
% 
% 
% ax_s2 = subplot(1,5,2,polaraxes);
% ph_2 = polarhistogram(theta_s2,36);
% 
% ax_s3 = subplot(1,5,3,polaraxes);
% ph_3 = polarhistogram(theta_s3,36);
% 
% ax_s4 = subplot(1,5,4,polaraxes);
% ph_4 = polarhistogram(theta_s4,36);
% 
% ax_s5 = subplot(1,5,5,polaraxes);
% ph_5 = polarhistogram(theta,36);
% 
% title(ax_s1,'Strategy 1');
% title(ax_s2,'Strategy 2');
% title(ax_s3,'Strategy 3');
% title(ax_s4,'Strategy 4');
% title(ax_s5,'Experiment');
% 
% rticklabels(ax_s1,{});
% rticklabels(ax_s2,{});
% rticklabels(ax_s3,{});
% rticklabels(ax_s4,{});
% rticklabels(ax_s5,{});
% ax_s1.RGrid = 'off';
% ax_s2.RGrid = 'off';
% ax_s3.RGrid = 'off';
% ax_s4.RGrid = 'off';
% ax_s5.RGrid = 'off';
% ax_s1.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s1.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ax_s2.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s2.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ax_s3.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s3.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ax_s4.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s4.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ax_s5.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s5.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ph_1.FaceColor = [0.2 0.2 0.2];
% ph_2.FaceColor = [0.2 0.2 0.2];
% ph_3.FaceColor = [0.2 0.2 0.2];
% ph_4.FaceColor = [0.2 0.2 0.2];
% ph_5.FaceColor = [0.2 0.2 0.2];

%% v=2
% fig = figure(2);
% fig.Color = 'w';
% fig.Position = [166.6000 281.8000 1.7008e+03 503.2000];
% ax_s1 = subplot(1,5,1,polaraxes);
% ph_1 = polarhistogram(theta_s1_02,36);
% 
% 
% ax_s2 = subplot(1,5,2,polaraxes);
% ph_2 = polarhistogram(theta_s2_02,36);
% 
% ax_s3 = subplot(1,5,3,polaraxes);
% ph_3 = polarhistogram(theta_s3_02,36);
% 
% ax_s4 = subplot(1,5,4,polaraxes);
% ph_4 = polarhistogram(theta_s4_02,36);
% 
% ax_s5 = subplot(1,5,5,polaraxes);
% ph_5 = polarhistogram(theta_02,36);
% 
% title(ax_s1,'Strategy 1');
% title(ax_s2,'Strategy 2');
% title(ax_s3,'Strategy 3');
% title(ax_s4,'Strategy 4');
% 
% rticklabels(ax_s1,{});
% rticklabels(ax_s2,{});
% rticklabels(ax_s3,{});
% rticklabels(ax_s4,{});
% rticklabels(ax_s5,{});
% ax_s1.RGrid = 'off';
% ax_s2.RGrid = 'off';
% ax_s3.RGrid = 'off';
% ax_s4.RGrid = 'off';
% ax_s5.RGrid = 'off';
% ax_s1.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s1.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ax_s2.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s2.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ax_s3.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s3.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ax_s4.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s4.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ax_s5.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s5.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ph_1.FaceColor = [0 0.2 0.7];
% ph_2.FaceColor = [0 0.2 0.7];
% ph_3.FaceColor = [0 0.2 0.7];
% ph_4.FaceColor = [0 0.2 0.7];
% ph_5.FaceColor = [0 0.2 0.7];

%% v=11
% fig = figure(3);
% fig.Color = 'w';
% fig.Position = [166.6000 281.8000 1.7008e+03 503.2000];
% ax_s1 = subplot(1,5,1,polaraxes);
% ph_1 = polarhistogram(theta_s1_11,36);
% 
% 
% ax_s2 = subplot(1,5,2,polaraxes);
% ph_2 = polarhistogram(theta_s2_11,36);
% 
% ax_s3 = subplot(1,5,3,polaraxes);
% ph_3 = polarhistogram(theta_s3_11,36);
% 
% ax_s4 = subplot(1,5,4,polaraxes);
% ph_4 = polarhistogram(theta_s4_11,36);
% 
% ax_s5 = subplot(1,5,5,polaraxes);
% ph_5 = polarhistogram(theta_11,36);
% 
% title(ax_s1,'Strategy 1');
% title(ax_s2,'Strategy 2');
% title(ax_s3,'Strategy 3');
% title(ax_s4,'Strategy 4');
% 
% rticklabels(ax_s1,{});
% rticklabels(ax_s2,{});
% rticklabels(ax_s3,{});
% rticklabels(ax_s4,{});
% rticklabels(ax_s5,{});
% ax_s1.RGrid = 'off';
% ax_s2.RGrid = 'off';
% ax_s3.RGrid = 'off';
% ax_s4.RGrid = 'off';
% ax_s5.RGrid = 'off';
% ax_s1.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s1.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ax_s2.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s2.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ax_s3.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s3.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ax_s4.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s4.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ax_s5.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s5.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ph_1.FaceColor = [0.7 0.2 0.2];
% ph_2.FaceColor = [0.7 0.2 0.2];
% ph_3.FaceColor = [0.7 0.2 0.2];
% ph_4.FaceColor = [0.7 0.2 0.2];
% ph_5.FaceColor = [0.7 0.2 0.2];

%% v=20
% fig = figure(4);
% fig.Color = 'w';
% fig.Position = [166.6000 281.8000 1.7008e+03 503.2000];
% ax_s1 = subplot(1,5,1,polaraxes);
% ph_1 = polarhistogram(theta_s1_20,36);
% 
% 
% ax_s2 = subplot(1,5,2,polaraxes);
% ph_2 = polarhistogram(theta_s2_20,36);
% 
% ax_s3 = subplot(1,5,3,polaraxes);
% ph_3 = polarhistogram(theta_s3_20,36);
% 
% ax_s4 = subplot(1,5,4,polaraxes);
% ph_4 = polarhistogram(theta_s4_20,36);
% 
% ax_s5 = subplot(1,5,5,polaraxes);
% ph_5 = polarhistogram(theta_20,36);
% 
% title(ax_s1,'Strategy 1');
% title(ax_s2,'Strategy 2');
% title(ax_s3,'Strategy 3');
% title(ax_s4,'Strategy 4');
% 
% 
% rticklabels(ax_s1,{});
% rticklabels(ax_s2,{});
% rticklabels(ax_s3,{});
% rticklabels(ax_s4,{});
% rticklabels(ax_s5,{});
% ax_s1.RGrid = 'off';
% ax_s2.RGrid = 'off';
% ax_s3.RGrid = 'off';
% ax_s4.RGrid = 'off';
% ax_s5.RGrid = 'off';
% ax_s1.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s1.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ax_s2.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s2.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ax_s3.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s3.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ax_s4.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s4.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ax_s5.ThetaTick = [0 45 90 135 180 225 270 315];
% ax_s5.ThetaTickLabel = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ph_1.FaceColor = [0.2 0.6 0.2];
% ph_2.FaceColor = [0.2 0.6 0.2];
% ph_3.FaceColor = [0.2 0.6 0.2];
% ph_4.FaceColor = [0.2 0.6 0.2];
% ph_5.FaceColor = [0.2 0.6 0.2];
%% histogram comparison between experiment and orthogonal/DO strategy per batch

figure('position',[276,534,1073,339]);c = get(gca,'colororder');
subplot(1,3,1);
histogram(abs(theta_OT_02),18,'BinLimits',[0,pi],'Normalization','Probability','facecolor',c(1,:));hold on;
histogram(abs(theta_02),18,'BinLimits',[0,pi],'Normalization','Probability','facecolor',c(5,:));
xlabel('turn angle $\theta$');
xticks([0,pi/2,pi]);xticklabels({'0','$\pi/2$','$\pi$'});
ylabel('probability');
legend('Orthogonal','experiment');
title('slow predator');
subplot(1,3,2);
histogram(abs(theta_OT_11),18,'BinLimits',[0,pi],'Normalization','Probability','facecolor',c(2,:));hold on;
histogram(abs(theta_11),18,'BinLimits',[0,pi],'Normalization','Probability','facecolor',c(5,:));
xlabel('turn angle $\theta$');
xticks([0,pi/2,pi]);xticklabels({'0','$\pi/2$','$\pi$'});
ylabel('probability');
legend('Orthogonal','experiment');
title('medium speed predator');
subplot(1,3,3);
histogram(abs(theta_OT_20),18,'BinLimits',[0,pi],'Normalization','Probability','facecolor',c(3,:));hold on;
histogram(abs(theta_20),18,'BinLimits',[0,pi],'Normalization','Probability','facecolor',c(5,:));
xlabel('turn angle $\theta$');
xticks([0,pi/2,pi]);xticklabels({'0','$\pi/2$','$\pi$'});
ylabel('probability');
legend('Orthogonal','experiment');
title('fast predator');

figure('position',[276,534,1073,339]);
subplot(1,3,1);
histogram(abs(theta_DOR_02),18,'BinLimits',[0,pi],'Normalization','Probability','facecolor',c(1,:));hold on;
histogram(abs(theta_02),18,'BinLimits',[0,pi],'Normalization','Probability','facecolor',c(5,:));
xlabel('turn angle $\theta$');
xticks([0,pi/2,pi]);xticklabels({'0','$\pi/2$','$\pi$'});
ylabel('probability');
legend('Weihs-Webb','experiment');
title('slow predator');
subplot(1,3,2);
histogram(abs(theta_DOR_11),18,'BinLimits',[0,pi],'Normalization','Probability','facecolor',c(2,:));hold on;
histogram(abs(theta_11),18,'BinLimits',[0,pi],'Normalization','Probability','facecolor',c(5,:));
xlabel('turn angle $\theta$');
xticks([0,pi/2,pi]);xticklabels({'0','$\pi/2$','$\pi$'});
ylabel('probability');
legend('Weihs-Webb','experiment');
title('medium speed predator');
subplot(1,3,3);
histogram(abs(theta_DOR_20),18,'BinLimits',[0,pi],'Normalization','Probability','facecolor',c(3,:));hold on;
histogram(abs(theta_20),18,'BinLimits',[0,pi],'Normalization','Probability','facecolor',c(5,:));
xlabel('turn angle $\theta$');
xticks([0,pi/2,pi]);xticklabels({'0','$\pi/2$','$\pi$'});
ylabel('probability');ylim([0,0.15])
legend('Weihs-Webb','experiment');
title('fast predator');
%%
load('turnVSEnergy.mat');
figure('position',[276,534,1073,339]);c = get(gca,'colororder');
subplot(1,3,1);
% yyaxis left
[N1,edges] = histcounts(abs(theta_OT_02),18,'BinLimits',[0,pi],'Normalization','Probability');
[N2,edges] = histcounts(abs(theta_02),18,'BinLimits',[0,pi],'Normalization','Probability');
% [N1_02,edges] = histcounts(abs(theta_s3_02),18,'BinLimits',[0,pi],'Normalization','Probability');
% [N2_02,edges] = histcounts(abs(theta_02),18,'BinLimits',[0,pi],'Normalization','Probability');
% [N1_11,edges] = histcounts(abs(theta_s3_11),18,'BinLimits',[0,pi],'Normalization','Probability');
% [N2_11,edges] = histcounts(abs(theta_11),18,'BinLimits',[0,pi],'Normalization','Probability');
% [N1_20,edges] = histcounts(abs(theta_s3_20),18,'BinLimits',[0,pi],'Normalization','Probability');
% [N2_20,edges] = histcounts(abs(theta_20),18,'BinLimits',[0,pi],'Normalization','Probability');
% hold on;
% xx = (edges(1:end-1)+edges(2:end))/2;
% yy1 = median([N1_02;N1_11;N1_20]);
% yy1_neg = yy1 - min([N1_02;N1_11;N1_20]);
% yy1_pos = max([N1_02;N1_11;N1_20]) - yy1;
% yy2 = median([N2_02;N2_11;N2_20]);
% yy2_neg = yy2 - min([N2_02;N2_11;N2_20]);
% yy2_pos = max([N2_02;N2_11;N2_20]) - yy2;
% colormap autumn
% C = [yy2 - yy1 flip(yy2)-flip(yy1)]*1';
% patch([xx flip(xx)]',[yy2 flip(yy1)]',C,'FaceColor','interp');
% e1=errorbar(xx,yy1,yy1_neg,yy1_pos,'o-','linewidth',1);
% e2=errorbar(xx,yy2,yy2_neg,yy2_pos,'o-','linewidth',1);
% colorbar
% xlabel('turn angle $\theta$');
% xlim([0,pi]);
% xticks([0,pi/2,pi]);xticklabels({'0','$\pi/2$','$\pi$'});
% ylabel('probability difference');
% legend([e1, e2],{'Orthogonal','experiment'});
% subplot(1,2,2);
% plot(turn,total_energySimple,'r-');
% ylabel('effort');
% xlabel('turn angle $\theta$');
% xlim([0,pi]);
% xticks([0,pi/2,pi]);xticklabels({'0','$\pi/2$','$\pi$'});

bar((edges(1:end-1)+edges(2:end))/2,N1-N2,'facecolor',c(1,:));
xlabel('turn angle $\theta$');
xticks([0,pi/2,pi]);xticklabels({'0','$\pi/2$','$\pi$'});
ylabel('error');
ylim([-0.065,0.065]);
% legend('Orthogonal','experiment');
yyaxis right
plot(turn,total_energySimple,'r-');
ylabel('effort');

title('slow predator');
subplot(1,3,2);
yyaxis left
[N1,edges] = histcounts(abs(theta_OT_11),18,'BinLimits',[0,pi],'Normalization','Probability');
[N2,edges] = histcounts(abs(theta_11),18,'BinLimits',[0,pi],'Normalization','Probability');
hold on;
bar((edges(1:end-1)+edges(2:end))/2,N1-N2,'facecolor',c(2,:));
xlabel('turn angle $\theta$');
xticks([0,pi/2,pi]);xticklabels({'0','$\pi/2$','$\pi$'});
ylabel('error');
ylim([-0.065,0.065]);
% legend('Orthogonal','experiment');
yyaxis right
plot(turn,total_energySimple,'r-');
ylabel('effort');
title('medium speed predator');
subplot(1,3,3);
yyaxis left
[N1,edges] = histcounts(abs(theta_OT_20),18,'BinLimits',[0,pi],'Normalization','Probability');
[N2,edges] = histcounts(abs(theta_20),18,'BinLimits',[0,pi],'Normalization','Probability');
hold on;
bar((edges(1:end-1)+edges(2:end))/2,N1-N2,'facecolor',c(3,:));
xlabel('turn angle $\theta$');
xticks([0,pi/2,pi]);xticklabels({'0','$\pi/2$','$\pi$'});
ylabel('error');
ylim([-0.065,0.065]);
% legend('Orthogonal','experiment');
title('fast predator');
yyaxis right
plot(turn,total_energySimple,'r-');
ylabel('effort');
%% all data compapre
load('turnVSEnergy.mat');
figure('position',[276 412 284 385]);c = get(gca,'colororder');
subplot(3,1,1);
histogram(abs(theta),18,'BinLimits',[0,pi],'facecolor',c(1,:));hold on;
xlabel('turn angle $\theta$');xlim([0-pi/18,pi+pi/18]);
xticks([0,pi/2,pi]);xticklabels({'$0$','$90$','$180$'});
ylabel('counts');
ylim([0,100]);
title('experiment');
subplot(3,1,2);
histogram(abs(theta_DO),18,'BinLimits',[0,pi],'facecolor',c(1,:));hold on;
xlabel('turn angle $\theta$');xlim([0-pi/18,pi+pi/18]);
xticks([0,pi/2,pi]);xticklabels({'$0$','$90$','$180$'});
ylabel('counts');
ylim([0,100]);
title('distance-optimal');
subplot(3,1,3);
histogram(abs(theta_OT),18,'BinLimits',[0,pi],'facecolor',c(1,:));hold on;
xlabel('turn angle $\theta$');xlim([0-pi/18,pi+pi/18]);
xticks([0,pi/2,pi]);xticklabels({'$0$','$90$','$180$'});
ylabel('counts');
ylim([0,100]);
title('orthogonal');

figure('position',[276,534,320,200]);c = get(gca,'colororder');
% subplot(2,1,1);
yyaxis left
[N1,edges] = histcounts(abs(theta_OT),18,'BinLimits',[0,pi]);
[N2,edges] = histcounts(abs(theta),18,'BinLimits',[0,pi]);
hold on;
% bar((edges(1:end-1)+edges(2:end))/2,N1-N2,'facecolor',c(2,:));
p1 = plot((edges(1:end-1)+edges(2:end))/2,N1-N2,'rs-','linewidth',0.5,'MarkerEdgeColor','none','MarkerFaceColor','r');
% m = '$\mathcal{D}$';
% text((edges(1:end-1)+edges(2:end))/2,N1-N2,m,'FontSize',10,'color','b');
xlabel('turn angle $\theta$');xlim([0-pi/18,pi+pi/18]);
xticks([0,pi/2,pi]);xticklabels({'0','$90$','$180$'});
ylabel('error');
% ylim([-0.05,0.05]);
% legend('Orthogonal','experiment');
yyaxis right
plot(turn,total_energySimple,'k-');
ylim([0,0.5]);
ylabel('effort');
% title('Orthogonal');
% subplot(2,1,2);
yyaxis left
[N1,edges] = histcounts(abs(theta_DO),18,'BinLimits',[0,pi]);
[N2,edges] = histcounts(abs(theta),18,'BinLimits',[0,pi]);
hold on;
% bar((edges(1:end-1)+edges(2:end))/2,N1-N2,'facecolor',c(3,:));
p2 = plot((edges(1:end-1)+edges(2:end))/2,N1-N2,'bo-','linewidth',0.5,'MarkerEdgeColor','none','MarkerFaceColor','b');
% m = '$\mathcal{O}$';
% text((edges(1:end-1)+edges(2:end))/2,N1-N2,m,'FontSize',10,'color','red');

xlabel('turn angle $\theta$');xlim([0-pi/18,pi+pi/18]);
xticks([0,pi/2,pi]);xticklabels({'0','$90$','$180$'});
ylabel('error');
% ylim([-0.05,0.05]);

% title('Distance-optimal');
% yyaxis right
% plot(turn,total_energySimple,'r-');
ylabel('effort');
%% response distribution
load('turnVSEnergy.mat');
figure('position',[276,534,1073,339]);
subplot(1,3,1);
yyaxis left
h = histogram(abs(theta_02(u_02<100)),18,'BinLimits',[0,pi],'Normalization','probability','facecolor',c(1,:));
hold on; xline(mean(theta_02(u_02<100)),'--k','linewidth',1.5);
xlabel('turn angle $\theta$');
ylabel('probability');ylim([0,0.25]);
yyaxis right
plot(turn,1./total_energySimple,'r-');
ylabel('1/effort');
title('slow predator');
subplot(1,3,2);
yyaxis left
h = histogram(abs(theta_11(u_11<100)),18,'BinLimits',[0,pi],'Normalization','probability','facecolor',c(2,:));
hold on; xline(mean(theta_11(u_11<100)),'--k','linewidth',1.5);
xlabel('turn angle $\theta$');ylabel('probability');ylim([0,0.25]);
yyaxis right
plot(turn,1./total_energySimple,'r-');
ylabel('1/effort');
title('medium speed predator');
subplot(1,3,3);
yyaxis left
h = histogram(abs(theta_20(u_20<100)),18,'BinLimits',[0,pi],'Normalization','probability','facecolor',c(3,:));
hold on; xline(mean(theta_20(u_20<100)),'--k','linewidth',1.5);
xlabel('turn angle $\theta$');
ylabel('probability');ylim([0,0.25]);
yyaxis right
plot(turn,1./total_energySimple,'r-');
ylabel('1/effort');
title('fast predator');

%%




% %% residues
% thetaRes_1 = wrapToPi(theta - theta_s1);
% thetaRes_2 = wrapToPi(theta - theta_s2);
% thetaRes_3 = wrapToPi(theta - theta_s3);
% thetaRes_4 = wrapToPi(theta - theta_s4);
% 
% figure,histogram(thetaRes_2,10);
% hold on;
% histogram(thetaRes_3,10);
% histogram(thetaRes_4,10);
% legend('2','3','4');

%% 3-D KL-Divergence per v
% k = 1;
% X = [theta_02; phi_02; lam_02]; % s1
% Y = [theta_s1_02; phi_02; lam_02]; % s1
% kl_s1(1) = KLDivergence(X,Y,k);
% X = [theta_11; phi_11; lam_11]; % s1
% Y = [theta_s1_11; phi_11; lam_11]; % s1
% kl_s1(2) = KLDivergence(X,Y,k);
% X = [theta_20; phi_20; lam_20 ]; % s1
% Y = [theta_s1_20; phi_20; lam_20]; % s1
% kl_s1(3) = KLDivergence(X,Y,k);
% 
% X = [theta_02; phi_02; lam_02]; % s2
% Y = [theta_s2_02; phi_02; lam_02]; % s2
% kl_s2(1) = KLDivergence(X,Y,k);
% X = [theta_11; phi_11; lam_11]; % s2
% Y = [theta_s2_11; phi_11; lam_11]; % s2
% kl_s2(2) = KLDivergence(X,Y,k);
% X = [theta_20; phi_20; lam_20]; % s2
% Y = [theta_s2_20; phi_20; lam_20]; % s2
% kl_s2(3) = KLDivergence(X,Y,k);
% 
% X = [theta_02; phi_02; lam_02]; % s3
% Y = [theta_s3_02; phi_02; lam_02]; % s3
% kl_s3(1) = KLDivergence(X,Y,k);
% X = [theta_11; phi_11; lam_11]; % s3
% Y = [theta_s3_11; phi_11; lam_11]; % s3
% kl_s3(2) = KLDivergence(X,Y,k);
% X = [theta_20; phi_20; lam_20]; % s3
% Y = [theta_s3_20; phi_20; lam_20]; % s3
% kl_s3(3) = KLDivergence(X,Y,k);
% 
% X = [theta_02; phi_02; lam_02]; % s4
% Y = [theta_s4_02; phi_02; lam_02]; % s4
% kl_s4(1) = KLDivergence(X,Y,k);
% X = [theta_11; phi_11; lam_11]; % s4
% Y = [theta_s4_11; phi_11; lam_11]; % s4
% kl_s4(2) = KLDivergence(X,Y,k);
% X = [theta_20; phi_20; lam_20]; % s4
% Y = [theta_s4_20; phi_20; lam_20]; % s4
% kl_s4(3) = KLDivergence(X,Y,k);
% 
% X = [theta; phi; lam]; % s1
% Y = [theta_s1; phi; lam]; % s1
% kl_s1(4) = KLDivergence(X,Y,k);
% Y = [theta_s2; phi; lam]; % s2
% kl_s2(4) = KLDivergence(X,Y,k);
% Y = [theta_s3; phi; lam]; % s3
% kl_s3(4) = KLDivergence(X,Y,k);
% Y = [theta_s4; phi; lam]; % s4
% kl_s4(4) = KLDivergence(X,Y,k);
% 
% K = [kl_s1;kl_s2;kl_s3;kl_s4];
% figure,bar(K);
% legend('v=2','v=11','v=20','all');
%% 4-D KL-Divergence Bootstrapping
% k = 1;
% 
% kl_his = zeros(15,4);
% for j = 1:200
% r = randperm(699,600);
% X = [theta(r); phi(r); lam(r); v(r)]; % s1
% Y = [theta_s1(r); phi(r); lam(r); v(r)]; % s1
% % X = theta(r); % s1
% % Y = theta_s1(r); % s1
% 
% kl_s1 = KLDivergence(X,Y,k);
% kl_his(j,1) = kl_s1;
% 
% r = randperm(699,600);
% X = [theta(r); phi(r); lam(r); v(r)];
% Y = [theta_s2(r); phi(r); lam(r); v(r)]; % s2
% % X = theta(r); % s2
% % Y = theta_s2(r); % s2
% kl_s2 = KLDivergence(X,Y,k);
% kl_his(j,2) = kl_s2;
% 
% r = randperm(699,600);
% X = [theta(r); phi(r); lam(r); v(r)];
% Y = [theta_s3(r); phi(r); lam(r); v(r)]; % s3
% % X = theta(r); % s3
% % Y = theta_s3(r); % s3
% kl_s3 = KLDivergence(X,Y,k);
% kl_his(j,3) = kl_s3;
% 
% r = randperm(699,600);
% X = [theta(r); phi(r); lam(r); v(r)];
% Y = [theta_s4(r); phi(r); lam(r); v(r)]; % s4
% % X = theta(r); % s4
% % Y = theta_s4(r); % s4
% 
% kl_s4 = KLDivergence(X,Y,k);
% kl_his(j,4) = kl_s4;
% end



%%
% figure,histogram(kl_his(:,1),10);
% hold on;
% histogram(kl_his(:,2),10);
% histogram(kl_his(:,3),10);
% histogram(kl_his(:,4),10);
% legend('1','2','3','4');
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
% %%
% [h1,p1,ci1,stats1] = ttest(kl_his(:,1));
% [h2,p2,ci2,stats2] = ttest(kl_his(:,2));
% [h3,p3,ci3,stats3] = ttest(kl_his(:,3));
% [h4,p4,ci4,stats4] = ttest(kl_his(:,4));
% %%
% [h23,p23,ci23,stats23] = ttest(kl_his(:,2),kl_his(:,3));
% [h24,p24,ci24,stats24] = ttest(kl_his(:,2),kl_his(:,4));
% [h34,p34,ci34,stats34] = ttest(kl_his(:,3),kl_his(:,4));

%% speed ratio inffered from Weihs-Webb strategy
chi_s4 = sign(lam).*(wrapToPi(phi + lam - theta + pi));
figure,histogram(chi_s4,20,'normalization','probability');
xlim([-pi,pi]);
xlabel('$\chi$ by inverse Distance-optimal');
ylabel('probability');
xticks([-pi,-pi/2,0,pi/2,pi]);
xticklabels({'$-\pi$','$-\pi/2$','0','$\pi/2$','$\pi$'});
%%
ratio = exp(-6:0.01:5);
chi = real(acos(ratio));
figure
semilogx(ratio,chi);
xlim([0,100]);
ylim([0,pi/2]);
xlabel('speed ratio $u/v$')
ylabel('$\chi = \arccos(u/v)$')
yticks([0,pi/4,pi/2]);
yticklabels({'0','$\pi/4$','$\pi/2$'});