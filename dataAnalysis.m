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

colorcode_strategy = [239,65,56;
    41,152,213;
    75,190,108;
    150,109,183;
    255,143,41]/255;
map = flip(gray);
load('evasionData.mat');
% u is pre-computed larvae average speed using u = \delta d/\delta t
% v is predator speed
% phi is predator angular position
% lambda is 'relative heading' of predator
% theta is the evasion direction
% d is the distance between predator and prey
phi = data.Phi(isfinite(data.U));
u   = data.U(isfinite(data.U));
v   = data.V(isfinite(data.U))*0.01; % original unit is cm/s, converted to m/s
lam = data.Lambda(isfinite(data.U));
theta = data.Delta(isfinite(data.U));
d = data.D(isfinite(data.U));

psii = wrapToPi(phi+lam+pi);
theta = wrapToPi(theta);
assert(all(lam),'some lambda is zero, which will make the sign function below be 0');
assert(all(phi-pi),'some phi equals to pi, which will make the sign function below be 0');
%% Contralateral
branch = sign(phi-pi);
theta_CL = branch*pi/2;
theta_CL= wrapToPi(theta_CL);
%% Antipodal
branch = sign(phi-pi);
theta_AP = phi + branch*pi;
theta_AP= wrapToPi(theta_AP);
%% Parallel
theta_PL = psii;
theta_PL= wrapToPi(theta_PL);
%% Distance-optimal, u/v = 0.5
branch = sign(lam);
chi = acos(0.5);
theta_DO = psii - branch*chi;

theta_DO= wrapToPi(theta_DO);
%% Orthogonal
branch = sign(lam);
chi = acos(0);
theta_OT = psii - branch*chi;

theta_OT= wrapToPi(theta_OT);
% %% Distance-optimal, u/v from data (real)
% branch = sign(lam);
% chi = real(acos(u./v));
% theta_DOR = phi + lam + pi - branch.*chi;
% theta_DOR= wrapTo2Pi(theta_DOR);
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
d_02 = d(v == 0.02);
d_11 = d(v == 0.11);
d_20 = d(v == 0.2);

psii_02 = psii(v == 0.02);
psii_11 = psii(v == 0.11);
psii_20 = psii(v == 0.2);
% theta_CL_02 = theta_CL(v == 0.02);
% theta_CL_11 = theta_CL(v == 0.11);
% theta_CL_20 = theta_CL(v == 0.2);
% theta_AP_02 = theta_AP(v == 0.02);
% theta_AP_11 = theta_AP(v == 0.11);
% theta_AP_20 = theta_AP(v == 0.2);
% theta_PL_02 = theta_PL(v == 0.02);
% theta_PL_11 = theta_PL(v == 0.11);
% theta_PL_20 = theta_PL(v == 0.2);
% theta_DO_02 = theta_DO(v == 0.02);
% theta_DO_11 = theta_DO(v == 0.11);
% theta_DO_20 = theta_DO(v == 0.2);
% theta_OT_02 = theta_OT(v == 0.02);
% theta_OT_11 = theta_OT(v == 0.11);
% theta_OT_20 = theta_OT(v == 0.2);
%% Correlations

[mnPhi,muPhi,mlPhi] = circ_mean(phi');
[mnPhi_02,muPhi_02,mlPhi_02] = circ_mean(phi_02');
[mnPhi_11,muPhi_11,mlPhi_11] = circ_mean(phi_11');
[mnPhi_20,muPhi_20,mlPhi_20] = circ_mean(phi_20');
[mnLam,muLam,mlLam] = circ_mean(lam');
[mnLam_02,muLam_02,mlLam_02] = circ_mean(lam_02');
[mnLam_11,muLam_11,mlLam_11] = circ_mean(lam_11');
[mnLam_20,muLam_20,mlLam_20] = circ_mean(lam_20');
[mnPsi,muPsi,mlPsi] = circ_mean(psii');
[mnPsi_02,muPsi_02,mlPsi_02] = circ_mean(psii_02');
[mnPsi_11,muPsi_11,mlPsi_11] = circ_mean(psii_11');
[mnPsi_20,muPsi_20,mlPsi_20] = circ_mean(psii_20');
[mnTheta,muTheta,mlTheta] = circ_mean(theta');
[mnTheta_02,muTheta_02,mlTheta_02] = circ_mean(theta_02');
[mnTheta_11,muTheta_11,mlTheta_11] = circ_mean(theta_11');
[mnTheta_20,muTheta_20,mlTheta_20] = circ_mean(theta_20');

stdPhi = circ_std(phi');
stdPhi_02 = circ_std(phi_02');
stdPhi_11 = circ_std(phi_11');
stdPhi_20 = circ_std(phi_20');
stdLam = circ_std(lam');
stdLam_02 = circ_std(lam_02');
stdLam_11 = circ_std(lam_11');
stdLam_20 = circ_std(lam_20');
stdPsi = circ_std(psii');
stdPsi_02 = circ_std(psii_02');
stdPsi_11 = circ_std(psii_11');
stdPsi_20 = circ_std(psii_20');
stdTheta = circ_std(theta');
stdTheta_02 = circ_std(theta_02');
stdTheta_11 = circ_std(theta_11');
stdTheta_20 = circ_std(theta_20');

% phi-lambda
[rhoPL, pvalPL] = circ_corrcc((phi)',(lam)');
[rhoPL_02, pvalPL_02] = circ_corrcc((phi_02)',(lam_02)');
[rhoPL_11, pvalPL_11] = circ_corrcc((phi_11)',(lam_11)');
[rhoPL_20, pvalPL_20] = circ_corrcc((phi_20)',(lam_20)');

% phi-theta
[rhoPT, pvalPT] = circ_corrcc(phi',theta');
[rhoPT_02, pvalPT_02] = circ_corrcc(phi_02',theta_02');
[rhoPT_11, pvalPT_11] = circ_corrcc(phi_11',theta_11');
[rhoPT_20, pvalPT_20] = circ_corrcc(phi_20',theta_20');
fitPT = circular_linearFitting(phi, theta);
fitPT_02 = circular_linearFitting(phi_02, theta_02);
fitPT_11 = circular_linearFitting(phi_11, theta_11);
fitPT_20 = circular_linearFitting(phi_20, theta_20);
% lambda-theta
[rhoLT, pvalLT] = circ_corrcc(lam',theta');
[rhoLT_02, pvalLT_02] = circ_corrcc(lam_02',theta_02');
[rhoLT_11, pvalLT_11] = circ_corrcc(lam_11',theta_11');
[rhoLT_20, pvalLT_20] = circ_corrcc(lam_20',theta_20');
fitLT = circular_linearFitting(lam, theta);
fitLT_02 = circular_linearFitting(lam_02(lam_02>0), theta_02(lam_02>0));
fitLT_11 = circular_linearFitting(lam_11, theta_11);
fitLT_20 = circular_linearFitting(lam_20(lam_20>0), theta_20(lam_20>0));
% psi-theta
[rhoPLT, pvalPLT] = circ_corrcc(psii',theta');
[rhoPLT_02, pvalPLT_02] = circ_corrcc(psii_02',theta_02');
[rhoPLT_11, pvalPLT_11] = circ_corrcc(psii_11',theta_11');
[rhoPLT_20, pvalPLT_20] = circ_corrcc(psii_20',theta_20');
fitPLT = circular_linearFitting(psii, theta);
fitPLT_02 = circular_linearFitting(psii_02, theta_02);
fitPLT_11 = circular_linearFitting(psii_11, theta_11);
fitPLT_20 = circular_linearFitting(psii_20, theta_20);
% psi-theta (lambda > 0)
[rhoPLT_pL, pvalPLT_pL] = circ_corrcc(psii(lam>0)',theta(lam>0)');
[rhoPLT_pL_02, pvalPLT_pL_02] = circ_corrcc(psii_02(lam_02>0)',theta_02(lam_02>0)');
[rhoPLT_pL_11, pvalPLT_pL_11] = circ_corrcc(psii_11(lam_11>0)',theta_11(lam_11>0)');
[rhoPLT_pL_20, pvalPLT_pL_20] = circ_corrcc(psii_20(lam_20>0)',theta_20(lam_20>0)');
fitPLT_pL = circular_linearFitting(psii(lam>0), theta(lam>0));
fitPLT_pL_02 = circular_linearFitting(psii_02(lam_02>0), theta_02(lam_02>0));
fitPLT_pL_11 = circular_linearFitting(psii_11(lam_11>0), theta_11(lam_11>0));
fitPLT_pL_20 = circular_linearFitting(psii_20(lam_20>0), theta_20(lam_20>0));
% psi-theta (lambda < 0)
[rhoPLT_nL, pvalPLT_nL] = circ_corrcc(psii(lam<0)',theta(lam<0)');
[rhoPLT_nL_02, pvalPLT_nL_02] = circ_corrcc(psii_02(lam_02<0)',theta_02(lam_02<0)');
[rhoPLT_nL_11, pvalPLT_nL_11] = circ_corrcc(psii_11(lam_11<0)',theta_11(lam_11<0)');
[rhoPLT_nL_20, pvalPLT_nL_20] = circ_corrcc(psii_20(lam_20<0)',theta_20(lam_20<0)');
fitPLT_nL = circular_linearFitting(psii(lam<0), theta(lam<0));
fitPLT_nL_02 = circular_linearFitting(psii_02(lam_02<0), theta_02(lam_02<0));
fitPLT_nL_11 = circular_linearFitting(psii_11(lam_11<0), theta_11(lam_11<0));
fitPLT_nL_20 = circular_linearFitting(psii_20(lam_20<0), theta_20(lam_20<0));
% d-theta
[rhoDT, pvalDT] = circ_corrcl(theta', d');
[rhoDT_02, pvalDT_02] = circ_corrcl(theta_02', d_02');
[rhoDT_11, pvalDT_11] = circ_corrcl(theta_11', d_11');
[rhoDT_20, pvalDT_20] = circ_corrcl(theta_20', d_20');

disp(['ϕ mean (with upper/lower 95% confidence interval):' newline num2str(mnPhi) '(' num2str(muPhi) ',' num2str(mlPhi) ')  '...
    num2str(mnPhi_02) '(' num2str(muPhi_02) ',' num2str(mlPhi_02) ')  '...
    num2str(mnPhi_11) '(' num2str(muPhi_11) ',' num2str(mlPhi_11) ')  '...
    num2str(mnPhi_20) '(' num2str(muPhi_20) ',' num2str(mlPhi_20) ')']);
disp(['ψ mean (with upper/lower 95% confidence interval):' newline num2str(mnPsi) '(' num2str(muPsi) ',' num2str(mlPsi) ')  '...
    num2str(mnPsi_02) '(' num2str(muPsi_02) ',' num2str(mlPsi_02) ')  '...
    num2str(mnPsi_11) '(' num2str(muPsi_11) ',' num2str(mlPsi_11) ')  '...
    num2str(mnPsi_20) '(' num2str(muPsi_20) ',' num2str(mlPsi_20) ')']);
disp(['λ mean (with upper/lower 95% confidence interval):' newline num2str(mnLam) '(' num2str(muLam) ',' num2str(mlLam) ')  '...
    num2str(mnLam_02) '(' num2str(muLam_02) ',' num2str(mlLam_02) ')  '...
    num2str(mnLam_11) '(' num2str(muLam_11) ',' num2str(mlLam_11) ')  '...
    num2str(mnLam_20) '(' num2str(muLam_20) ',' num2str(mlLam_20) ')']);
disp(['θ mean (with upper/lower 95% confidence interval):' newline num2str(mnTheta) '(' num2str(muTheta) ',' num2str(mlTheta) ')  '...
    num2str(mnTheta_02) '(' num2str(muTheta_02) ',' num2str(mlTheta_02) ')  '...
    num2str(mnTheta_11) '(' num2str(muTheta_11) ',' num2str(mlTheta_11) ')  '...
    num2str(mnTheta_20) '(' num2str(muTheta_20) ',' num2str(mlTheta_20) ')']);

disp(['d mean & std:' newline num2str(mean(d_02)) ' & ' num2str(std(d_02)) ' & ' num2str(mean(d_11)) ' & ' num2str(std(d_11)) ' & ' num2str(mean(d_20)) ' & ' num2str(std(d_20)) ' & ' num2str(mean(d)) ' & ' num2str(std(d))]);
disp(['ϕ mean & std:' newline num2str(wrapTo2Pi(mnPhi_02*180/pi)) ' & ' num2str(stdPhi_02*180/pi) ' & ' num2str(wrapTo2Pi(mnPhi_11*180/pi)) ' & ' num2str(stdPhi_11*180/pi) ' & ' num2str(wrapTo2Pi(mnPhi_20*180/pi)) ' & ' num2str(stdPhi_20*180/pi) ' & ' num2str(wrapTo2Pi(mnPhi*180/pi)) ' & ' num2str(stdPhi*180/pi)]);
disp(['ψ mean & std:' newline num2str(wrapTo2Pi(mnPsi_02*180/pi)) ' & ' num2str(stdPsi_02*180/pi) ' & ' num2str(wrapTo2Pi(mnPsi_11*180/pi)) ' & ' num2str(stdPsi_11*180/pi) ' & ' num2str(wrapTo2Pi(mnPsi_20*180/pi)) ' & ' num2str(stdPsi_20*180/pi) ' & ' num2str(wrapTo2Pi(mnPsi*180/pi)) ' & ' num2str(stdPsi*180/pi)]);
disp(['λ mean & std:' newline num2str(wrapTo2Pi(mnLam_02*180/pi)) ' & ' num2str(stdLam_02*180/pi) ' & ' num2str(wrapTo2Pi(mnLam_11*180/pi)) ' & ' num2str(stdLam_11*180/pi) ' & ' num2str(wrapTo2Pi(mnLam_20*180/pi)) ' & ' num2str(stdLam_20*180/pi) ' & ' num2str(wrapTo2Pi(mnLam*180/pi)) ' & ' num2str(stdLam*180/pi)]);
disp(['θ mean & std:' newline num2str(wrapTo2Pi(mnTheta_02*180/pi)) ' & ' num2str(stdTheta_02*180/pi) ' & ' num2str(wrapTo2Pi(mnTheta_11*180/pi)) ' & ' num2str(stdTheta_11*180/pi) ' & ' num2str(wrapTo2Pi(mnTheta_20*180/pi)) ' & ' num2str(stdTheta_20*180/pi) ' & ' num2str(wrapTo2Pi(mnTheta*180/pi)) ' & ' num2str(stdTheta*180/pi)]);

disp(['d & θ correlation & regression:' newline ...
    num2str(rhoDT_02) ' & ' num2str(pvalDT_02) ' & ' ' ' ' & ' num2str(rhoDT_11) ' & ' num2str(pvalDT_11) ' & ' ' ' ' & ' num2str(rhoDT_20) ' & ' num2str(pvalDT_20) ' & ' ' ' ' & ' num2str(rhoDT) ' & ' num2str(pvalDT) ' & ' ' ']);
disp(['ϕ & θ correlation & regression:' newline ...
    num2str(rhoPT_02) ' & ' num2str(pvalPT_02) ' & ' num2str(fitPT_02(1)) ' & ' num2str(rhoPT_11) ' & ' num2str(pvalPT_11) ' & ' num2str(fitPT_11(1)) ' & ' num2str(rhoPT_20) ' & ' num2str(pvalPT_20) ' & ' num2str(fitPT_20(1)) ' & ' num2str(rhoPT) ' & ' num2str(pvalPT) ' & ' num2str(fitPT(1))]);
disp(['ψ & θ correlation & regression:' newline ...
    num2str(rhoPLT_02) ' & ' num2str(pvalPLT_02) ' & ' num2str(fitPLT_02(1)) ' & ' num2str(rhoPLT_11) ' & ' num2str(pvalPLT_11) ' & ' num2str(fitPLT_11(1)) ' & ' num2str(rhoPLT_20) ' & ' num2str(pvalPLT_20) ' & ' num2str(fitPLT_20(1)) ' & ' num2str(rhoPLT) ' & ' num2str(pvalPLT) ' & ' num2str(fitPLT(1))]);
disp(['λ & θ correlation & regression:' newline ...
    num2str(rhoLT_02) ' & ' num2str(pvalLT_02) ' & ' num2str(fitLT_02(1)) ' & ' num2str(rhoLT_11) ' & ' num2str(pvalLT_11) ' & ' num2str(fitLT_11(1)) ' & ' num2str(rhoLT_20) ' & ' num2str(pvalLT_20) ' & ' num2str(fitLT_20(1)) ' & ' num2str(rhoLT) ' & ' num2str(pvalLT) ' & ' num2str(fitLT(1))]);
disp(['ψ (λ>0) & θ correlation & regression:' newline ...
    num2str(rhoPLT_pL_02) ' & ' num2str(pvalPLT_pL_02) ' & ' num2str(fitPLT_pL_02(1)) ' & ' num2str(rhoPLT_pL_11) ' & ' num2str(pvalPLT_pL_11) ' & ' num2str(fitPLT_pL_11(1)) ' & ' num2str(rhoPLT_pL_20) ' & ' num2str(pvalPLT_pL_20) ' & ' num2str(fitPLT_pL_20(1)) ' & ' num2str(rhoPLT_pL) ' & ' num2str(pvalPLT_pL) ' & ' num2str(fitPLT_pL(1))]);
disp(['ψ (λ<0) & θ correlation & regression:' newline ...
    num2str(rhoPLT_nL_02) ' & ' num2str(pvalPLT_nL_02) ' & ' num2str(fitPLT_nL_02(1)) ' & ' num2str(rhoPLT_nL_11) ' & ' num2str(pvalPLT_nL_11) ' & ' num2str(fitPLT_nL_11(1)) ' & ' num2str(rhoPLT_nL_20) ' & ' num2str(pvalPLT_nL_20) ' & ' num2str(fitPLT_nL_20(1)) ' & ' num2str(rhoPLT_nL) ' & ' num2str(pvalPLT_nL) ' & ' num2str(fitPLT_nL(1))]);
disp('phi-lam');
disp(['correlation: ' num2str([rhoPL,rhoPL_02,rhoPL_11,rhoPL_20])]);
disp(['p-value: ' num2str([pvalPL,pvalPL_02,pvalPL_11,pvalPL_20])]);

disp('phi-theta');
disp(['correlation: ' num2str([rhoPT,rhoPT_02,rhoPT_11,rhoPT_20])]);
disp(['p-value: ' num2str([pvalPT,pvalPT_02,pvalPT_11,pvalPT_20])]);

disp('lam-theta');
disp(['correlation: ' num2str([rhoLT,rhoLT_02,rhoLT_11,rhoLT_20])]);
disp(['p-value: ' num2str([pvalLT,pvalLT_02,pvalLT_11,pvalLT_20])]);

disp('phi+lam+pi,psi-theta');
disp(['correlation: ' num2str([rhoPLT,rhoPLT_02,rhoPLT_11,rhoPLT_20])]);
disp(['p-value: ' num2str([pvalPLT,pvalPLT_02,pvalPLT_11,pvalPLT_20])]);

disp('phi+lam+pi,psi-theta(lam>0)');
disp(['correlation: ' num2str([rhoPLT_pL,rhoPLT_pL_02,rhoPLT_pL_11,rhoPLT_pL_20])]);
disp(['p-value: ' num2str([pvalPLT_pL,pvalPLT_pL_02,pvalPLT_pL_11,pvalPLT_pL_20])]);

disp('phi+lam+pi,psi-theta(lam<0)');
disp(['correlation: ' num2str([rhoPLT_nL,rhoPLT_nL_02,rhoPLT_nL_11,rhoPLT_nL_20])]);
disp(['p-value: ' num2str([pvalPLT_nL,pvalPLT_nL_02,pvalPLT_nL_11,pvalPLT_nL_20])]);

disp('d-theta');
disp(['correlation: ' num2str([rhoDT,rhoDT_02,rhoDT_11,rhoDT_20])]);
disp(['p-value: ' num2str([pvalDT,pvalDT_02,pvalDT_11,pvalDT_20])]);
%% Start from here, using pre-run results
load('correlation_analysis_results_temp.mat');
%% phi-lambda
f = figure('Position',[361 676 1140 437]);
colormap(map);
xedges = 0:pi/15:2*pi; yedges = -pi:pi/15:pi;
[X,Y] = meshgrid(xedges,yedges);
hc_exp = histcounts2(wrapTo2Pi(phi_02),wrapToPi(lam_02),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,1), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi$');ylabel('$\lambda$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});


hc_exp = histcounts2(wrapTo2Pi(phi_11),wrapToPi(lam_11),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,2), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi$');ylabel('$\lambda$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(wrapTo2Pi(phi_20),wrapToPi(lam_20),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,3), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi$');ylabel('$\lambda$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(wrapTo2Pi(phi),wrapToPi(lam),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,4), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi$'); ylabel('$\lambda$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

%% phi-theta
f = figure('Position',[361 676 1140 437]);
colormap(map);
xedges = 0:pi/15:2*pi; yedges = -pi:pi/15:pi;
[X,Y] = meshgrid(xedges,yedges);
hc_exp = histcounts2(wrapTo2Pi(phi_02),theta_02,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,1), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(phi_02), max(phi_02)];
hold on; plot(x_ends, x_ends*fitPT_02(1) + fitPT_02(2),'r');
axis equal;
xlabel('$\phi$');ylabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});


hc_exp = histcounts2(wrapTo2Pi(phi_11),theta_11,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,2), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(phi_11), max(phi_11)];
hold on; plot(x_ends, x_ends*fitPT_11(1) + fitPT_11(2),'r');
axis equal;
xlabel('$\phi$');ylabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(wrapTo2Pi(phi_20),theta_20,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,3), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(phi_20), max(phi_20)];
hold on; plot(x_ends, x_ends*fitPT_20(1) + fitPT_20(2),'r');
axis equal;
xlabel('$\phi$');ylabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(wrapTo2Pi(phi),theta,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,4), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(phi), max(phi)];
hold on; plot(x_ends, x_ends*fitPT(1) + fitPT(2),'r');
axis equal;
xlabel('$\phi$');ylabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

%% lambda-theta
f = figure('Position',[361 676 1140 437]);
colormap(map);
xedges = -pi:pi/15:pi; yedges = -pi:pi/15:pi;
[X,Y] = meshgrid(xedges,yedges);
hc_exp = histcounts2(wrapToPi(lam_02),theta_02,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,1), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(lam_02), max(lam_02)];
hold on; plot(x_ends, x_ends*fitLT_02(1) + fitLT_02(2),'r');
axis equal;
xlabel('$\lambda$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([-pi,0,pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});


hc_exp = histcounts2(wrapToPi(lam_11),theta_11,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,2), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(lam_11), max(lam_11)];
hold on; plot(x_ends, x_ends*fitLT_11(1) + fitLT_11(2),'r');
axis equal;
xlabel('$\lambda$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([-pi,0,pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(wrapToPi(lam_20),theta_20,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,3), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(lam_20), max(lam_20)];
hold on; plot(x_ends, x_ends*fitLT_20(1) + fitLT_20(2),'r');
axis equal;
xlabel('$\lambda$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([-pi,0,pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(wrapToPi(lam),theta,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,4), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(lam), max(lam)];
hold on; plot(x_ends, x_ends*fitLT(1) + fitLT(2),'r');
axis equal;
xlabel('$\lambda$'); ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([-pi,0,pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

%% phi+lambda (psi)-theta
f = figure('Position',[361 676 1140 437]);
colormap(map);
xedges = -pi:pi/15:pi; yedges = -pi:pi/15:pi;
[X,Y] = meshgrid(xedges,yedges);
hc_exp = histcounts2(psii_02,theta_02,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,1), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(psii_02), max(psii_02)];
hold on; plot(x_ends, x_ends*fitPLT_02(1) + fitPLT_02(2),'r');
axis equal;
xlabel('$\psi$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([-pi,0,pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});



hc_exp = histcounts2(psii_11,theta_11,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,2), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(psii_11), max(psii_11)];
hold on; plot(x_ends, x_ends*fitPLT_11(1) + fitPLT_11(2),'r');
axis equal;
xlabel('$\psi$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([-pi,0,pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(psii_20,theta_20,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,3), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(psii_20), max(psii_20)];
hold on; plot(x_ends, x_ends*fitPLT_20(1) + fitPLT_20(2),'r');
axis equal;
xlabel('$\psi$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([-pi,0,pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(psii,theta,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,4), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(psii), max(psii)];
hold on; plot(x_ends, x_ends*fitPLT(1) + fitPLT(2),'r');
axis equal;
xlabel('$\psi$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([-pi,0,pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});
%% phi - psi - theta
f = figure('Position',[361 676 500 500]);
colormap(map);
markerSize = 12;
subplot(1,3,1), sc = scatter3(phi_02, psii_02, theta_02,markerSize);
axis equal;
xlabel('$\psi$');ylabel('$\psi$');zlabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);zlim([-pi,pi])
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);zticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});
yticklabels({'$-\pi$','$0$','$\pi$'});
zticklabels({'$-\pi$','$0$','$\pi$'});


subplot(1,3,2), sc = scatter3(phi_11, psii_11, theta_11,markerSize);
axis equal;
xlabel('$\psi$');ylabel('$\psi$');zlabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);zlim([-pi,pi])
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);zticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});
yticklabels({'$-\pi$','$0$','$\pi$'});
zticklabels({'$-\pi$','$0$','$\pi$'});
subplot(1,3,3), sc = scatter3(phi_20, psii_20, theta_20,markerSize);
axis equal;
xlabel('$\psi$');ylabel('$\psi$');zlabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);zlim([-pi,pi])
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);zticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});
yticklabels({'$-\pi$','$0$','$\pi$'});
zticklabels({'$-\pi$','$0$','$\pi$'});
%% |phi|-|lambda|
f = figure('Position',[361 676 1140 437]);
colormap(map);
xedges = 0:pi/30:pi; yedges = 0:pi/30:pi;
[X,Y] = meshgrid(xedges,yedges);
hc_exp = histcounts2(abs(wrapToPi(phi_02)),abs(lam_02),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,3,1), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi$');ylabel('$\lambda$');
xlim([0,pi]);ylim([0,pi]);
% axis off; % colorbar;
xticks([0,pi]);yticks([0,pi]);
xticklabels({'$0$','$\pi$'});yticklabels({'$0$','$\pi$'});

hc_exp = histcounts2(abs(wrapToPi(phi_11)),abs(lam_11),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,3,2), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi$');ylabel('$\lambda$');
xlim([0,pi]);ylim([0,pi]);
% axis off; % colorbar;
xticks([0,pi]);yticks([0,pi]);
xticklabels({'$0$','$\pi$'});yticklabels({'$0$','$\pi$'});

hc_exp = histcounts2(abs(wrapToPi(phi_20)),abs(lam_20),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,3,3), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi$');ylabel('$\lambda$');
xlim([0,pi]);ylim([0,pi]);
% axis off; % colorbar;
xticks([0,pi]);yticks([0,pi]);
xticklabels({'$0$','$\pi$'});yticklabels({'$0$','$\pi$'});

%% phi+lambda (psi)-theta (lambda>0/lambda<0)
f = figure('Position',[361 676 1140 800]);
colormap(map);
xedges = -pi:pi/15:pi; yedges = -pi:pi/15:pi;
[X,Y] = meshgrid(xedges,yedges);
hc_exp = histcounts2(psii_02(lam_02>0),theta_02(lam_02>0),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,4,1), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(psii_02(lam_02>0)), max(psii_02(lam_02>0))];
hold on; plot(x_ends, x_ends*fitPLT_pL_02(1) + fitPLT_pL_02(2),'r');
axis equal;
xlabel('$\psi$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(psii_11(lam_11>0),theta_11(lam_11>0),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,4,2), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(psii_11(lam_11>0)), max(psii_11(lam_11>0))];
hold on; plot(x_ends, x_ends*fitPLT_pL_11(1) + fitPLT_pL_11(2),'r');
axis equal;
xlabel('$\psi$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([-pi,0,pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(psii_20(lam_20>0),theta_20(lam_20>0),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,4,3), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(psii_20(lam_20>0)), max(psii_20(lam_20>0))];
hold on; plot(x_ends, x_ends*fitPLT_pL_20(1) + fitPLT_pL_20(2),'r');
axis equal;
xlabel('$\psi$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([-pi,0,pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(psii(lam>0),theta(lam>0),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,4,4), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(psii(lam>0)), max(psii(lam>0))];
hold on; plot(x_ends, x_ends*fitPLT_pL(1) + fitPLT_pL(2),'r');
axis equal;
xlabel('$\psi$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([-pi,0,pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});


hc_exp = histcounts2(psii_02(lam_02<0),theta_02(lam_02<0),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,4,5), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(psii_02(lam_02<0)), max(psii_02(lam_02<0))];
hold on; plot(x_ends, x_ends*fitPLT_nL_02(1) + fitPLT_nL_02(2),'r');
axis equal;
xlabel('$\psi$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
xticks([-pi,0,pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(psii_11(lam_11<0),theta_11(lam_11<0),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,4,6), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(psii_11(lam_11<0)), max(psii_11(lam_11<0))];
hold on; plot(x_ends, x_ends*fitPLT_nL_11(1) + fitPLT_nL_11(2),'r');
axis equal;
xlabel('$\psi$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(psii_20(lam_20<0),theta_20(lam_20<0),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,4,7), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(psii_20(lam_20<0)), max(psii_20(lam_20<0))];
hold on; plot(x_ends, x_ends*fitPLT_nL_20(1) + fitPLT_nL_20(2),'r');
axis equal;
xlabel('$\psi$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
xticks([-pi,0,pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(psii(lam<0),theta(lam<0),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,4,8), h_exp = pcolor(X,Y,hc_exp');
x_ends = [min(psii(lam<0)), max(psii(lam<0))];
hold on; plot(x_ends, x_ends*fitPLT_nL_20(1) + fitPLT_nL_20(2),'r');
axis equal;
xlabel('$\psi$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
xticks([-pi,0,pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});
%% d-theta
f = figure('Position',[361 676 500 437]);
colormap(map);
scale = 2*pi/0.06;
xedges = (0:0.002:0.06)*scale; yedges = -pi:pi/15:pi;
[X,Y] = meshgrid(xedges,yedges);
hc_exp = histcounts2(d_02*scale, theta_02,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,1), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$d$');ylabel('$\theta$');
xlim([0, 0.06]*scale);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,0.03,0.06]*scale);yticks([-pi,0,pi]);
xticklabels({'$0$','$3$','$6$'});yticklabels({'$-\pi$','$0$','$\pi$'});


hc_exp = histcounts2(d_11*scale, theta_11,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,2), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$d$');ylabel('$\theta$');
xlim([0, 0.06]*scale);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,0.03,0.06]*scale);yticks([-pi,0,pi]);
xticklabels({'$0$','$3$','$6$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(d_20*scale,theta_20,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,3), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$d$');ylabel('$\theta$');
xlim([0, 0.06]*scale);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,0.03,0.06]*scale);yticks([-pi,0,pi]);
xticklabels({'$0$','$3$','$6$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(d*scale,theta,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,4,4), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$d$');ylabel('$\theta$');
xlim([0, 0.06]*scale);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,0.03,0.06]*scale);yticks([-pi,0,pi]);
xticklabels({'$0$','$3$','$6$'});yticklabels({'$-\pi$','$0$','$\pi$'});
%% phi, no noise

fp = figure();fp.Position = [338 321 500 180];
colormap(map);
% map = zeros(300,3);
% map(1:230,2) = linspace(1,0,230);
% map(1:230,3) = linspace(1,0,230);
% map(1:230,1) = 1;
% map(231:300,1) =linspace(1,0.65,70);
% colormap(map);

[X,Y] = meshgrid(0:pi/15:2*pi,-pi:pi/15:pi);
xedges = 0:pi/15:2*pi; yedges = -pi:pi/15:pi;

hc_exp = histcounts2(phi,theta,xedges,yedges,'normalization','probability');
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,4,1), h_exp = pcolor(X,Y,hc_exp');
axis equal;title('$\phi$,no noise');axis off;
xticks([0,2*pi]);yticks([-pi,pi]);
colorbar;
hc_CL = histcounts2(phi,theta_CL,xedges,yedges,'normalization','probability');
hc_CL = [hc_CL zeros(size(hc_CL,1),1)];hc_CL = [hc_CL;zeros(1,size(hc_CL,2))];
subplot(2,4,2), h_s1 = pcolor(X,Y,hc_CL');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);
colorbar;
hc_AP = histcounts2(phi,theta_AP,xedges,yedges,'normalization','probability');
hc_AP = [hc_AP zeros(size(hc_AP,1),1)];hc_AP = [hc_AP;zeros(1,size(hc_AP,2))];
subplot(2,4,3), h_s2 = pcolor(X,Y,hc_AP');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);
colorbar;
hc_PL = histcounts2(phi,theta_PL,xedges,yedges,'normalization','probability');
hc_PL = [hc_PL zeros(size(hc_PL,1),1)];hc_PL = [hc_PL;zeros(1,size(hc_PL,2))];
subplot(2,4,4), h_s3 = pcolor(X,Y,hc_PL');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);
colorbar;
hc_DO = histcounts2(phi,theta_DO,xedges,yedges,'normalization','probability');
hc_DO = [hc_DO zeros(size(hc_DO,1),1)];hc_DO = [hc_DO;zeros(1,size(hc_DO,2))];
subplot(2,4,6), h_s4_1 = pcolor(X,Y,hc_DO');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);
colorbar;
hc_OT = histcounts2(phi,theta_OT,xedges,yedges,'normalization','probability');
hc_OT = [hc_OT zeros(size(hc_OT,1),1)];hc_OT = [hc_OT;zeros(1,size(hc_OT,2))];
subplot(2,4,7), h_s4_2 = pcolor(X,Y,hc_OT');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);
colorbar;

%% lambda, no noise
fl = figure();fl.Position = [338 321 500 180];
colormap(map);
[X,Y] = meshgrid(-pi:pi/15:pi,-pi:pi/15:pi);
xedges = -pi:pi/15:pi; yedges = -pi:pi/15:pi;
hc_exp = histcounts2(lam,theta,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,4,1), h_exp = pcolor(X,Y,hc_exp');title('$\lambda$,no noise');
axis equal; % colorbar;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_CL = histcounts2(lam,theta_CL,xedges,yedges);
hc_CL = [hc_CL zeros(size(hc_CL,1),1)];hc_CL = [hc_CL;zeros(1,size(hc_CL,2))];
subplot(2,4,2), h_s1 = pcolor(X,Y,hc_CL');
axis equal;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_AP = histcounts2(lam,theta_AP,xedges,yedges);
hc_AP = [hc_AP zeros(size(hc_AP,1),1)];hc_AP = [hc_AP;zeros(1,size(hc_AP,2))];
subplot(2,4,3), h_s2 = pcolor(X,Y,hc_AP');
axis equal;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_PL = histcounts2(lam,theta_PL,xedges,yedges);
hc_PL = [hc_PL zeros(size(hc_PL,1),1)];hc_PL = [hc_PL;zeros(1,size(hc_PL,2))];
subplot(2,4,4), h_s3 = pcolor(X,Y,hc_PL');
axis equal;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_DO = histcounts2(lam,theta_DO,xedges,yedges);
hc_DO = [hc_DO zeros(size(hc_DO,1),1)];hc_DO = [hc_DO;zeros(1,size(hc_DO,2))];
subplot(2,4,6), h_s4_1 = pcolor(X,Y,hc_DO');
axis equal;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_OT = histcounts2(lam,theta_OT,xedges,yedges);
hc_OT = [hc_OT zeros(size(hc_OT,1),1)];hc_OT = [hc_OT;zeros(1,size(hc_OT,2))];
subplot(2,4,7), h_s4_2 = pcolor(X,Y,hc_OT');
axis equal;
xticks([-pi,pi]);yticks([-pi,pi]);

%% ADDING NOISE
load('opt_initialize.mat')
%% Contralateral

sig_phi = eta0_CL(1);
sig_theta = eta0_CL(2);

N = length(theta);
theta_CL_noise = zeros(1,N);
for i = 1:N

    phih = randraw('vonmises', [phi(i),1/sig_phi^2], [1 1]);
    thetah = sign(phih-pi)*pi/2;
    theta_CL_noise(i) = wrapToPi(randraw('vonmises', [thetah,1/sig_theta^2], [1 1]));

end
%% Antipodal

sig_phi = eta0_AP(1);
sig_theta = eta0_AP(1);

N = length(theta);
theta_AP_noise = zeros(1,N);
for i = 1:N

    phih = randraw('vonmises', [phi(i),1/sig_phi^2], [1 1]);
    thetah = phih + sign(pi-phih)*pi;
    theta_AP_noise(i) = wrapToPi(randraw('vonmises', [thetah,1/sig_theta^2], [1 1]));

end
%% Parallel
chi = 0;
sig_phi = eta0_PL(1);
sig_lam = eta0_PL(1);
sig_theta = eta0_PL(1);

theta_PL_noise = zeros(1,N);

for i = 1:N
    phih = randraw('vonmises', [phi(i),1/sig_phi^2], [1,1]);
    lamh = randraw('vonmises', [lam(i),1/sig_lam^2], [1,1]);
    thetah = phih + lamh + pi - sign(lamh)*chi;
    theta_PL_noise(i) = wrapToPi(randraw('vonmises', [thetah,1/sig_theta^2], [1,1]));
end

%% Distance-optimal
chi = acos(0.5);
sig_phi = eta0_DO(1);
sig_lam = eta0_DO(2);
sig_theta = eta0_DO(1);

theta_DO_noise = zeros(1,N);

for i = 1:N
    phih = randraw('vonmises', [phi(i),1/sig_phi^2], [1,1]);
    lamh = randraw('vonmises', [lam(i),1/sig_lam^2], [1,1]);
    thetah = phih + lamh + pi - sign(lamh)*chi;
    theta_DO_noise(i) = wrapToPi(randraw('vonmises', [thetah,1/sig_theta^2], [1,1]));
end

%% Orthogonal
chi = acos(0);
sig_phi = eta0_OT(1);
sig_lam = eta0_OT(2);
sig_theta = eta0_OT(1);

theta_OT_noise = zeros(1,N);

for i = 1:N
    phih = randraw('vonmises', [phi(i),1/sig_phi^2], [1,1]);
    lamh = randraw('vonmises', [lam(i),1/sig_lam^2], [1,1]);
    thetah = phih + lamh + pi - sign(lamh)*chi;
    theta_OT_noise(i) = wrapToPi(randraw('vonmises', [thetah,1/sig_theta^2], [1,1]));
end
%% phi, noise
fpn = figure();fpn.Position = [338 321 500 180];
colormap(map);
[X,Y] = meshgrid(0:pi/15:2*pi,-pi:pi/15:pi);
xedges = 0:pi/15:2*pi; yedges = -pi:pi/15:pi;
hc_exp = histcounts2(phi,theta,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,4,1), h_exp = pcolor(X,Y,hc_exp');
axis equal;axis off;title('$\phi$,noise');
xticks([0,2*pi]);yticks([-pi,pi]);

hc_CL = histcounts2(phi,theta_CL_noise,xedges,yedges);
hc_CL = [hc_CL zeros(size(hc_CL,1),1)];hc_CL = [hc_CL;zeros(1,size(hc_CL,2))];
subplot(2,4,2), h_s1 = pcolor(X,Y,hc_CL');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);

hc_AP = histcounts2(phi,theta_AP_noise,xedges,yedges);
hc_AP = [hc_AP zeros(size(hc_AP,1),1)];hc_AP = [hc_AP;zeros(1,size(hc_AP,2))];
subplot(2,4,3), h_s2 = pcolor(X,Y,hc_AP');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);


hc_DO = histcounts2(phi,theta_PL_noise,xedges,yedges);
hc_DO = [hc_DO zeros(size(hc_DO,1),1)];hc_DO = [hc_DO;zeros(1,size(hc_DO,2))];
subplot(2,4,6), h_s4_1 = pcolor(X,Y,hc_DO');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);

hc_OT = histcounts2(phi,theta_DO_noise,xedges,yedges);
hc_OT = [hc_OT zeros(size(hc_OT,1),1)];hc_OT = [hc_OT;zeros(1,size(hc_OT,2))];
subplot(2,4,7), h_s4_2 = pcolor(X,Y,hc_OT');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);

hc_s4_3 = histcounts2(phi,theta_OT_noise,xedges,yedges);
hc_s4_3 = [hc_s4_3 zeros(size(hc_s4_3,1),1)];hc_s4_3 = [hc_s4_3;zeros(1,size(hc_s4_3,2))];
subplot(2,4,8), h_s4_3 = pcolor(X,Y,hc_s4_3');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);



%% lambda, noise
fln = figure();fln.Position = [338 321 500 180];
colormap(map);
[X,Y] = meshgrid(0:pi/15:2*pi,-pi:pi/15:pi);
xedges = -pi:pi/15:pi; yedges = -pi:pi/15:pi;
hc_exp = histcounts2(lam,theta,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,4,1), h_exp = pcolor(X,Y,hc_exp');
axis equal;axis off;title('$\lambda$,noise');
xticks([-pi,pi]);yticks([-pi,pi]);

hc_CL = histcounts2(lam,theta_CL_noise,xedges,yedges);
hc_CL = [hc_CL zeros(size(hc_CL,1),1)];hc_CL = [hc_CL;zeros(1,size(hc_CL,2))];
subplot(2,4,2), h_s1 = pcolor(X,Y,hc_CL');
axis equal;axis off;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_AP = histcounts2(lam,theta_AP_noise,xedges,yedges);
hc_AP = [hc_AP zeros(size(hc_AP,1),1)];hc_AP = [hc_AP;zeros(1,size(hc_AP,2))];
subplot(2,4,3), h_s2 = pcolor(X,Y,hc_AP');
axis equal;axis off;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_DO = histcounts2(lam,theta_PL_noise,xedges,yedges);
hc_DO = [hc_DO zeros(size(hc_DO,1),1)];hc_DO = [hc_DO;zeros(1,size(hc_DO,2))];
subplot(2,4,6), h_s4_1 = pcolor(X,Y,hc_DO');
axis equal;axis off;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_OT = histcounts2(lam,theta_DO_noise,xedges,yedges);
hc_OT = [hc_OT zeros(size(hc_OT,1),1)];hc_OT = [hc_OT;zeros(1,size(hc_OT,2))];
subplot(2,4,7), h_s4_2 = pcolor(X,Y,hc_OT');
axis equal;axis off;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_s4_3 = histcounts2(lam,theta_OT_noise,xedges,yedges);
hc_s4_3 = [hc_s4_3 zeros(size(hc_s4_3,1),1)];hc_s4_3 = [hc_s4_3;zeros(1,size(hc_s4_3,2))];
subplot(2,4,8), h_s4_3 = pcolor(X,Y,hc_s4_3');
axis equal;axis off;
xticks([-pi,pi]);yticks([-pi,pi]);
% figure, polarhistogram(theta_s3,36,'Normalization','pdf');
% subplot(2,4,5), h_s4 = histogram2(phi,theta,[20 20],'DisplayStyle','tile','ShowEmptyBins','on');
% axis equal;h_exp.XBinLimits = [0,2*pi];h_exp.YBinLimits = [-pi,pi];
% xticks([0,2*pi]);yticks([-pi,pi]);