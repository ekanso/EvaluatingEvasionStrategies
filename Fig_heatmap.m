clear;
close all;

load('evasionData.mat');
set(groot,'defaultLineLineWidth',2);
set(groot,'defaultFigureColor','w');
set(groot,'defaultTextFontsize',10);
set(groot,'defaultAxesFontsize',10);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesLineWidth',1);
phi = data.Phi(isfinite(data.U));
u   = data.U(isfinite(data.U));
v   = data.V(isfinite(data.U))*0.01;
lam = data.Lambda(isfinite(data.U));
theta = data.Delta(isfinite(data.U));

theta = wrapToPi(theta);

map = flip(gray);
%% S1
theta_s1 = sign(phi-pi)*pi/2;

theta_s1= wrapToPi(theta_s1);

%% S2
theta_s2 = phi + sign(pi-phi)*pi;

theta_s2= wrapToPi(theta_s2);
%% S3
theta_s3 = phi + lam + pi - ((lam>=0)-0.5)*pi;

theta_s3= wrapToPi(theta_s3);
%% S4,u/v=1
chi = 0;
theta_s4_1 = phi + lam + pi - ((lam>=0)-0.5)*chi*2;

theta_s4_1= wrapToPi(theta_s4_1);
%% S4,u/v=0.2
chi = acos(0.2);
theta_s4_2 = phi + lam + pi - ((lam>=0)-0.5)*chi*2;

theta_s4_2= wrapToPi(theta_s4_2);
%% S4,u/v=0.1
chi = acos(0.1);
theta_s4_3 = phi + lam + pi - ((lam>=0)-0.5)*chi*2;

theta_s4_3= wrapToPi(theta_s4_3);
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

theta_s1_02 = theta_s1(v == 0.02);
theta_s1_11 = theta_s1(v == 0.11);
theta_s1_20 = theta_s1(v == 0.2);
theta_s2_02 = theta_s2(v == 0.02);
theta_s2_11 = theta_s2(v == 0.11);
theta_s2_20 = theta_s2(v == 0.2);
theta_s3_02 = theta_s3(v == 0.02);
theta_s3_11 = theta_s3(v == 0.11);
theta_s3_20 = theta_s3(v == 0.2);
theta_s4_1_02 = theta_s4_1(v == 0.02);
theta_s4_1_11 = theta_s4_1(v == 0.11);
theta_s4_1_20 = theta_s4_1(v == 0.2);
theta_s4_2_02 = theta_s4_2(v == 0.02);
theta_s4_2_11 = theta_s4_2(v == 0.11);
theta_s4_2_20 = theta_s4_2(v == 0.2);
theta_s4_3_02 = theta_s4_3(v == 0.02);
theta_s4_3_11 = theta_s4_3(v == 0.11);
theta_s4_3_20 = theta_s4_3(v == 0.2);
%% phi-lambda
f = figure();
% fl.Position = [338 321 500 180];
colormap(map);
xedges = 0:pi/15:2*pi; yedges = -pi:pi/15:pi;
[X,Y] = meshgrid(xedges,yedges);
hc_exp = histcounts2(wrapTo2Pi(phi_02),wrapToPi(lam_02),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,3,1), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi$');ylabel('$\lambda$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});


hc_exp = histcounts2(wrapTo2Pi(phi_11),wrapToPi(lam_11),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,3,2), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi$');ylabel('$\lambda$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(wrapTo2Pi(phi_20),wrapToPi(lam_20),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,3,3), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi$');ylabel('$\lambda$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

%% phi-theta
f = figure();
% fl.Position = [338 321 500 180];
colormap(map);
xedges = 0:pi/15:2*pi; yedges = -pi:pi/15:pi;
[X,Y] = meshgrid(xedges,yedges);
hc_exp = histcounts2(wrapTo2Pi(phi_02),theta_02,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,3,1), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi$');ylabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});


hc_exp = histcounts2(wrapTo2Pi(phi_11),theta_11,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,3,2), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi$');ylabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(wrapTo2Pi(phi_20),theta_20,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,3,3), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi$');ylabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

%% lambda-theta
f = figure();
% fl.Position = [338 321 500 180];
colormap(map);
xedges = -pi:pi/15:pi; yedges = -pi:pi/15:pi;
[X,Y] = meshgrid(xedges,yedges);
hc_exp = histcounts2(wrapToPi(lam_02),theta_02,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,3,1), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\lambda$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([-pi,0,pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});


hc_exp = histcounts2(wrapToPi(lam_11),theta_11,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,3,2), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\lambda$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([-pi,0,pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(wrapToPi(lam_20),theta_20,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,3,3), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\lambda$');ylabel('$\theta$');
xlim([-pi,pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([-pi,0,pi]);yticks([-pi,0,pi]);
xticklabels({'$-\pi$','$0$','$\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

%% phi+lambda-theta
f = figure();
% fl.Position = [338 321 500 180];
colormap(map);
xedges = 0:pi/15:2*pi; yedges = -pi:pi/15:pi;
[X,Y] = meshgrid(xedges,yedges);
hc_exp = histcounts2(wrapTo2Pi(phi_02+lam_02),theta_02,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,3,1), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi+\lambda$');ylabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});



hc_exp = histcounts2(wrapTo2Pi(phi_11+lam_11),theta_11,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,3,2), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi+\lambda$');ylabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(wrapTo2Pi(phi_20+lam_20),theta_20,xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(1,3,3), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi+\lambda$');ylabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});
%% |phi|-|lambda|
f = figure();
% fl.Position = [338 321 500 180];
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

%% phi+lambda-theta (lambda>0/lambda<0)
f = figure();
% fl.Position = [338 321 500 180];
colormap(map);
xedges = 0:pi/15:2*pi; yedges = -pi:pi/15:pi;
[X,Y] = meshgrid(xedges,yedges);
hc_exp = histcounts2(wrapTo2Pi(phi_02(lam_02>0)+lam_02(lam_02>0)),theta_02(lam_02>0),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,3,1), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi+\lambda$');ylabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(wrapTo2Pi(phi_02(lam_02<0)+lam_02(lam_02<0)),theta_02(lam_02<0),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,3,4), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi+\lambda$');ylabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(wrapTo2Pi(phi_11(lam_11>0)+lam_11(lam_11>0)),theta_11(lam_11>0),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,3,2), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi+\lambda$');ylabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(wrapTo2Pi(phi_11(lam_11<0)+lam_11(lam_11<0)),theta_11(lam_11<0),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,3,5), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi+\lambda$');ylabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(wrapTo2Pi(phi_20(lam_20>0)+lam_20(lam_20>0)),theta_20(lam_20>0),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,3,3), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi+\lambda$');ylabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);
% axis off; % colorbar;
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});

hc_exp = histcounts2(wrapTo2Pi(phi_20(lam_20<0)+lam_20(lam_20<0)),theta_20(lam_20<0),xedges,yedges);
hc_exp = [hc_exp zeros(size(hc_exp,1),1)];hc_exp = [hc_exp;zeros(1,size(hc_exp,2))];
subplot(2,3,6), h_exp = pcolor(X,Y,hc_exp');
axis equal;
xlabel('$\phi+\lambda$');ylabel('$\theta$');
xlim([0,2*pi]);ylim([-pi,pi]);
xticks([0,pi,2*pi]);yticks([-pi,0,pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});yticklabels({'$-\pi$','$0$','$\pi$'});
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
hc_s1 = histcounts2(phi,theta_s1,xedges,yedges,'normalization','probability');
hc_s1 = [hc_s1 zeros(size(hc_s1,1),1)];hc_s1 = [hc_s1;zeros(1,size(hc_s1,2))];
subplot(2,4,2), h_s1 = pcolor(X,Y,hc_s1');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);
colorbar;
hc_s2 = histcounts2(phi,theta_s2,xedges,yedges,'normalization','probability');
hc_s2 = [hc_s2 zeros(size(hc_s2,1),1)];hc_s2 = [hc_s2;zeros(1,size(hc_s2,2))];
subplot(2,4,3), h_s2 = pcolor(X,Y,hc_s2');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);
colorbar;
hc_s3 = histcounts2(phi,theta_s3,xedges,yedges,'normalization','probability');
hc_s3 = [hc_s3 zeros(size(hc_s3,1),1)];hc_s3 = [hc_s3;zeros(1,size(hc_s3,2))];
subplot(2,4,4), h_s3 = pcolor(X,Y,hc_s3');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);
colorbar;
hc_s4_1 = histcounts2(phi,theta_s4_1,xedges,yedges,'normalization','probability');
hc_s4_1 = [hc_s4_1 zeros(size(hc_s4_1,1),1)];hc_s4_1 = [hc_s4_1;zeros(1,size(hc_s4_1,2))];
subplot(2,4,6), h_s4_1 = pcolor(X,Y,hc_s4_1');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);
colorbar;
hc_s4_2 = histcounts2(phi,theta_s4_2,xedges,yedges,'normalization','probability');
hc_s4_2 = [hc_s4_2 zeros(size(hc_s4_2,1),1)];hc_s4_2 = [hc_s4_2;zeros(1,size(hc_s4_2,2))];
subplot(2,4,7), h_s4_2 = pcolor(X,Y,hc_s4_2');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);
colorbar;
hc_s4_3 = histcounts2(phi,theta_s4_3,xedges,yedges,'normalization','probability');
hc_s4_3 = [hc_s4_3 zeros(size(hc_s4_3,1),1)];hc_s4_3 = [hc_s4_3;zeros(1,size(hc_s4_3,2))];
subplot(2,4,8), h_s4_3 = pcolor(X,Y,hc_s4_3');
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

hc_s1 = histcounts2(lam,theta_s1,xedges,yedges);
hc_s1 = [hc_s1 zeros(size(hc_s1,1),1)];hc_s1 = [hc_s1;zeros(1,size(hc_s1,2))];
subplot(2,4,2), h_s1 = pcolor(X,Y,hc_s1');
axis equal;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_s2 = histcounts2(lam,theta_s2,xedges,yedges);
hc_s2 = [hc_s2 zeros(size(hc_s2,1),1)];hc_s2 = [hc_s2;zeros(1,size(hc_s2,2))];
subplot(2,4,3), h_s2 = pcolor(X,Y,hc_s2');
axis equal;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_s3 = histcounts2(lam,theta_s3,xedges,yedges);
hc_s3 = [hc_s3 zeros(size(hc_s3,1),1)];hc_s3 = [hc_s3;zeros(1,size(hc_s3,2))];
subplot(2,4,4), h_s3 = pcolor(X,Y,hc_s3');
axis equal;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_s4_1 = histcounts2(lam,theta_s4_1,xedges,yedges);
hc_s4_1 = [hc_s4_1 zeros(size(hc_s4_1,1),1)];hc_s4_1 = [hc_s4_1;zeros(1,size(hc_s4_1,2))];
subplot(2,4,6), h_s4_1 = pcolor(X,Y,hc_s4_1');
axis equal;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_s4_2 = histcounts2(lam,theta_s4_2,xedges,yedges);
hc_s4_2 = [hc_s4_2 zeros(size(hc_s4_2,1),1)];hc_s4_2 = [hc_s4_2;zeros(1,size(hc_s4_2,2))];
subplot(2,4,7), h_s4_2 = pcolor(X,Y,hc_s4_2');
axis equal;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_s4_3 = histcounts2(lam,theta_s4_3,xedges,yedges);
hc_s4_3 = [hc_s4_3 zeros(size(hc_s4_3,1),1)];hc_s4_3 = [hc_s4_3;zeros(1,size(hc_s4_3,2))];
subplot(2,4,8), h_s4_3 = pcolor(X,Y,hc_s4_3');
axis equal;
xticks([-pi,pi]);yticks([-pi,pi]);
% ADDING NOISE
load('optNLLs.mat')
%% S1

sig_phi = etaCL(1);
sig_theta = etaCL(2);

N = length(theta);
theta_s1_n = zeros(1,N);
for i = 1:N

    phih = randraw('vonmises', [phi(i),1/sig_phi^2], [1 1]);
    thetah = sign(phih-pi)*pi/2;
    theta_s1_n(i) = wrapToPi(randraw('vonmises', [thetah,1/sig_theta^2], [1 1]));

end
%% S2

sig_phi = etaAP(1);
sig_theta = etaAP(2);

N = length(theta);
theta_s2_n = zeros(1,N);
for i = 1:N

    phih = randraw('vonmises', [phi(i),1/sig_phi^2], [1 1]);
    thetah = phih + sign(pi-phih)*pi;
    theta_s2_n(i) = wrapToPi(randraw('vonmises', [thetah,1/sig_theta^2], [1 1]));

end
%% Parallel
chi = 0;
sig_phi = etaPL(1);
sig_lam = etaPL(2);
sig_theta = etaPL(3);

theta_s4_1_n = zeros(1,N);

for i = 1:N
    phih = randraw('vonmises', [phi(i),1/sig_phi^2], [1,1]);
    lamh = randraw('vonmises', [lam(i),1/sig_lam^2], [1,1]);
    thetah = phih + lamh + pi - sign(lamh)*chi;
    theta_s4_1_n(i) = wrapToPi(randraw('vonmises', [thetah,1/sig_theta^2], [1,1]));
end

%% Distance-optimal
chi = acos(0.5);
sig_phi = etaDO(1);
sig_lam = etaDO(2);
sig_theta = etaDO(3);

theta_s4_2_n = zeros(1,N);

for i = 1:N
    phih = randraw('vonmises', [phi(i),1/sig_phi^2], [1,1]);
    lamh = randraw('vonmises', [lam(i),1/sig_lam^2], [1,1]);
    thetah = phih + lamh + pi - sign(lamh)*chi;
    theta_s4_2_n(i) = wrapToPi(randraw('vonmises', [thetah,1/sig_theta^2], [1,1]));
end

%% Orthogonal
chi = acos(0);
sig_phi = etaOT(1);
sig_lam = etaOT(2);
sig_theta = etaOT(3);

theta_s4_3_n = zeros(1,N);

for i = 1:N
    phih = randraw('vonmises', [phi(i),1/sig_phi^2], [1,1]);
    lamh = randraw('vonmises', [lam(i),1/sig_lam^2], [1,1]);
    thetah = phih + lamh + pi - sign(lamh)*chi;
    theta_s4_3_n(i) = wrapToPi(randraw('vonmises', [thetah,1/sig_theta^2], [1,1]));
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

hc_s1 = histcounts2(phi,theta_s1_n,xedges,yedges);
hc_s1 = [hc_s1 zeros(size(hc_s1,1),1)];hc_s1 = [hc_s1;zeros(1,size(hc_s1,2))];
subplot(2,4,2), h_s1 = pcolor(X,Y,hc_s1');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);

hc_s2 = histcounts2(phi,theta_s2_n,xedges,yedges);
hc_s2 = [hc_s2 zeros(size(hc_s2,1),1)];hc_s2 = [hc_s2;zeros(1,size(hc_s2,2))];
subplot(2,4,3), h_s2 = pcolor(X,Y,hc_s2');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);


hc_s4_1 = histcounts2(phi,theta_s4_1_n,xedges,yedges);
hc_s4_1 = [hc_s4_1 zeros(size(hc_s4_1,1),1)];hc_s4_1 = [hc_s4_1;zeros(1,size(hc_s4_1,2))];
subplot(2,4,6), h_s4_1 = pcolor(X,Y,hc_s4_1');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);

hc_s4_2 = histcounts2(phi,theta_s4_2_n,xedges,yedges);
hc_s4_2 = [hc_s4_2 zeros(size(hc_s4_2,1),1)];hc_s4_2 = [hc_s4_2;zeros(1,size(hc_s4_2,2))];
subplot(2,4,7), h_s4_2 = pcolor(X,Y,hc_s4_2');
axis equal;axis off;
xticks([0,2*pi]);yticks([-pi,pi]);

hc_s4_3 = histcounts2(phi,theta_s4_3_n,xedges,yedges);
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

hc_s1 = histcounts2(lam,theta_s1_n,xedges,yedges);
hc_s1 = [hc_s1 zeros(size(hc_s1,1),1)];hc_s1 = [hc_s1;zeros(1,size(hc_s1,2))];
subplot(2,4,2), h_s1 = pcolor(X,Y,hc_s1');
axis equal;axis off;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_s2 = histcounts2(lam,theta_s2_n,xedges,yedges);
hc_s2 = [hc_s2 zeros(size(hc_s2,1),1)];hc_s2 = [hc_s2;zeros(1,size(hc_s2,2))];
subplot(2,4,3), h_s2 = pcolor(X,Y,hc_s2');
axis equal;axis off;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_s4_1 = histcounts2(lam,theta_s4_1_n,xedges,yedges);
hc_s4_1 = [hc_s4_1 zeros(size(hc_s4_1,1),1)];hc_s4_1 = [hc_s4_1;zeros(1,size(hc_s4_1,2))];
subplot(2,4,6), h_s4_1 = pcolor(X,Y,hc_s4_1');
axis equal;axis off;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_s4_2 = histcounts2(lam,theta_s4_2_n,xedges,yedges);
hc_s4_2 = [hc_s4_2 zeros(size(hc_s4_2,1),1)];hc_s4_2 = [hc_s4_2;zeros(1,size(hc_s4_2,2))];
subplot(2,4,7), h_s4_2 = pcolor(X,Y,hc_s4_2');
axis equal;axis off;
xticks([-pi,pi]);yticks([-pi,pi]);

hc_s4_3 = histcounts2(lam,theta_s4_3_n,xedges,yedges);
hc_s4_3 = [hc_s4_3 zeros(size(hc_s4_3,1),1)];hc_s4_3 = [hc_s4_3;zeros(1,size(hc_s4_3,2))];
subplot(2,4,8), h_s4_3 = pcolor(X,Y,hc_s4_3');
axis equal;axis off;
xticks([-pi,pi]);yticks([-pi,pi]);
% figure, polarhistogram(theta_s3,36,'Normalization','pdf');
% subplot(2,4,5), h_s4 = histogram2(phi,theta,[20 20],'DisplayStyle','tile','ShowEmptyBins','on');
% axis equal;h_exp.XBinLimits = [0,2*pi];h_exp.YBinLimits = [-pi,pi];
% xticks([0,2*pi]);yticks([-pi,pi]);