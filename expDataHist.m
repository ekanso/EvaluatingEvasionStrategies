clear;
close all;
set(groot,'defaultLineLineWidth',2);
set(groot,'defaultFigureColor','w');
set(groot,'defaultTextFontsize',7);
set(groot,'defaultAxesFontsize',7);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesLineWidth',1);
load('evasionData.mat');
%%
phi = data.Phi(isfinite(data.U));
u   = data.U(isfinite(data.U));
v   = data.V(isfinite(data.U))*0.01;
lam = data.Lambda(isfinite(data.U));
theta = data.Delta(isfinite(data.U));
d = data.D(isfinite(data.U));

phi = wrapTo2Pi(phi);
lam = wrapToPi(lam);
theta = wrapToPi(theta);
psii = wrapToPi(phi+lam+pi);


phi_02 = phi(v == 0.02);
phi_11 = phi(v == 0.11);
phi_20 = phi(v == 0.2);
lam_02 = lam(v == 0.02);
lam_11 = lam(v == 0.11);
lam_20 = lam(v == 0.2);
psii_02 = psii(v == 0.02);
psii_11 = psii(v == 0.11);
psii_20 = psii(v == 0.2);
u_02 = u(v == 0.02);
u_11 = u(v == 0.11);
u_20 = u(v == 0.2);
d_02 = d(v == 0.02);
d_11 = d(v == 0.11);
d_20 = d(v == 0.2);
theta_02 = theta(v == 0.02);
theta_11 = theta(v == 0.11);
theta_20 = theta(v == 0.2);


%%
% in the order of (slow, mid-speed, fast, all) datasets
f = figure;
subplot(5,4,1),histogram(d_02,'Orientation','Vertical','Numbins',18,'BinLimits',[0,0.06]);

xticks([0,0.03,0.06]);
xlim([0,0.06]);
ylim([0,45]);
xlabel('$d$');
hold on, plot(mean(d_02)*[1,1],ylim,'-b')
subplot(5,4,2),histogram(d_11,'Orientation','Vertical','Numbins',18,'BinLimits',[0,0.06]);
xticks([0,0.03,0.06]);
xlim([0,0.06]);
ylim([0,45]);
xlabel('$d$');
hold on, plot(mean(d_11)*[1,1],ylim,'-b')
subplot(5,4,3),histogram(d_20,'Orientation','Vertical','Numbins',18,'BinLimits',[0,0.06]);
xticks([0,0.03,0.06]);
xlim([0,0.06]);
ylim([0,45]);
xlabel('$d$');
hold on, plot(mean(d_20)*[1,1],ylim,'-b')
subplot(5,4,4),histogram(d,'Orientation','Vertical','Numbins',18,'BinLimits',[0,0.06]);
xticks([0,0.03,0.06]);
xlim([0,0.06]);
ylim([0,100]);
xlabel('$d$');
hold on, plot(mean(d)*[1,1],ylim,'-b')

subplot(5,4,5),histogram(phi_02,'Orientation','Vertical','Numbins',18,'BinLimits',[0,2*pi]);
xticks([0,pi,2*pi]);
xlim([0,2*pi]);
ylim([0,40]);
xlabel('$\phi$');
hold on, plot(wrapTo2Pi(circ_mean(phi_02'))*[1,1],ylim,'-b')
subplot(5,4,6),histogram(phi_11,'Orientation','Vertical','Numbins',18,'BinLimits',[0,2*pi]);
xticks([0,pi,2*pi]);
xlim([0,2*pi]);
ylim([0,40]);
xlabel('$\phi$');
hold on, plot(wrapTo2Pi(circ_mean(phi_11'))*[1,1],ylim,'-b')
subplot(5,4,7),histogram(phi_20,'Orientation','Vertical','Numbins',18,'BinLimits',[0,2*pi]);
xticks([0,pi,2*pi]);
xlim([0,2*pi]);
ylim([0,40]);
xlabel('$\phi$');
hold on, plot(wrapTo2Pi(circ_mean(phi_20'))*[1,1],ylim,'-b')
subplot(5,4,8),histogram(phi,'Orientation','Vertical','Numbins',18,'BinLimits',[0,2*pi]);
xticks([0,pi,2*pi]);
xlim([0,2*pi]);
ylim([0,100]);
xlabel('$\phi$');
hold on, plot(wrapTo2Pi(circ_mean(phi'))*[1,1],ylim,'-b')

subplot(5,4,9),histogram(psii_02,'Orientation','Vertical','Numbins',18,'BinLimits',[-pi,pi]);
xticks([-pi,0,pi]);
xlim([-pi,pi]);
ylim([0,40]);
xlabel('$\psi$');
hold on, plot(wrapToPi(circ_mean(psii_02'))*[1,1],ylim,'-b')
subplot(5,4,10),histogram(psii_11,'Orientation','Vertical','Numbins',18,'BinLimits',[-pi,pi]);
xticks([-pi,0,pi]);
xlim([-pi,pi]);
ylim([0,40]);
xlabel('$\psi$');
hold on, plot(wrapToPi(circ_mean(psii_11'))*[1,1],ylim,'-b')
subplot(5,4,11),histogram(psii_20,'Orientation','Vertical','Numbins',18,'BinLimits',[-pi,pi]);
xticks([-pi,0,pi]);
xlim([-pi,pi]);
ylim([0,40]);
xlabel('$\psi$');
hold on, plot(wrapToPi(circ_mean(psii_20'))*[1,1],ylim,'-b')
subplot(5,4,12),histogram(psii,'Orientation','Vertical','Numbins',18,'BinLimits',[-pi,pi]);
xticks([-pi,0,pi]);
xlim([-pi,pi]);
ylim([0,100]);
xlabel('$\psi$');
hold on, plot(wrapToPi(circ_mean(psii'))*[1,1],ylim,'-b')

subplot(5,4,13),histogram(lam_02,'Orientation','Vertical','Numbins',18,'BinLimits',[-pi,pi]);
xticks([-pi,0,pi]);
xlim([-pi,pi]);
ylim([0,40]);
xlabel('$\lambda$');
hold on, plot(wrapToPi(circ_mean(lam_02'))*[1,1],ylim,'-b')
subplot(5,4,14),histogram(lam_11,'Orientation','Vertical','Numbins',18,'BinLimits',[-pi,pi]);
xticks([-pi,0,pi]);
xlim([-pi,pi]);
ylim([0,40]);
xlabel('$\lambda$');
hold on, plot(wrapToPi(circ_mean(lam_11'))*[1,1],ylim,'-b')
subplot(5,4,15),histogram(lam_20,'Orientation','Vertical','Numbins',18,'BinLimits',[-pi,pi]);
xticks([-pi,0,pi]);
xlim([-pi,pi]);
ylim([0,40]);
xlabel('$\lambda$');
hold on, plot(wrapToPi(circ_mean(lam_20'))*[1,1],ylim,'-b')
subplot(5,4,16),histogram(lam,'Orientation','Vertical','Numbins',18,'BinLimits',[-pi,pi]);
xticks([-pi,0,pi]);
xlim([-pi,pi]);
ylim([0,100]);
xlabel('$\lambda$');
hold on, plot(wrapToPi(circ_mean(lam'))*[1,1],ylim,'-b')


subplot(5,4,17),histogram(theta_02,'Orientation','Vertical','Numbins',18,'BinLimits',[-pi,pi]);
xticks([-pi,0,pi]);
xlim([-pi,pi]);
ylim([0,40]);
xlabel('$\theta$');
hold on, plot(wrapToPi(circ_mean(theta_02'))*[1,1],ylim,'-b')
subplot(5,4,18),histogram(theta_11,'Orientation','Vertical','Numbins',18,'BinLimits',[-pi,pi]);
xticks([-pi,0,pi]);
xlim([-pi,pi]);
ylim([0,40]);
xlabel('$\theta$');
hold on, plot(wrapToPi(circ_mean(theta_11'))*[1,1],ylim,'-b')
subplot(5,4,19),histogram(theta_20,'Orientation','Vertical','Numbins',18,'BinLimits',[-pi,pi]);
xticks([-pi,0,pi]);
xlim([-pi,pi]);
ylim([0,40]);
xlabel('$\theta$');
hold on, plot(wrapToPi(circ_mean(theta_20'))*[1,1],ylim,'-b')
subplot(5,4,20),histogram(theta,'Orientation','Vertical','Numbins',18,'BinLimits',[-pi,pi]);
xticks([-pi,0,pi]);
xlim([-pi,pi]);
ylim([0,100]);
xlabel('$\theta$');
hold on, plot(wrapToPi(circ_mean(theta'))*[1,1],ylim,'-b')

f.Position = [435 3 800 400];