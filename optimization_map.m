clear
close all
clc
set(groot,'defaultLineLineWidth',2);
set(groot,'defaultFigureColor','w');
set(groot,'defaultTextFontsize',24);
set(groot,'defaultAxesFontsize',24);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesLineWidth',1);
% global indx;
% indx = 1;
%% Load Data

load('evasionData.mat');
% load('s3_opt_result.mat');

%% Sort Data

phi = data.Phi(isfinite(data.U));
u   = data.U(isfinite(data.U));
v   = data.V(isfinite(data.U))*0.01;
lam = data.Lambda(isfinite(data.U));
theta = data.Delta(isfinite(data.U));

phi_02 = phi(v == 0.02);
phi_11 = phi(v == 0.11);
phi_20 = phi(v == 0.2);
lam_02 = lam(v == 0.02);
lam_11 = lam(v == 0.11);
lam_20 = lam(v == 0.2);
u_02 = u(v == 0.02);
u_11 = u(v == 0.11);
u_20 = u(v == 0.2);
v_02 = v(v == 0.02);
v_11 = v(v == 0.11);
v_20 = v(v == 0.2);
theta_02 = theta(v == 0.02);
theta_11 = theta(v == 0.11);
theta_20 = theta(v == 0.2);


%% Setup Data

Ns_02 = length(phi_02);

x_02 = [phi_02;...
    lam_02;...
    u_02;...
    v_02;...
    theta_02];
Ns_11 = length(phi_11);

x_11 = [phi_11;...
    lam_11;...
    u_11;...
    v_11;...
    theta_11];
Ns_20 = length(phi_20);

x_20 = [phi_20;...
    lam_20;...
    u_20;...
    v_20;...
    theta_20];

%% Setup Data

Ns = length(phi);

x = [phi;
    lam;
    u;
    v;
    theta];
%%
N = 120;
sig_theta = linspace(0.01, pi-0.01, N);
sig_phi = linspace(0.01, pi-0.01, N);
sig_lam = linspace(0.01, pi-0.01, N);
map = zeros(120);

for i = 1:N
    tic
    for j = 1:N
        eta = [sig_phi(i),sig_theta(j)];
        map(i,j) = NLL_CL_reimagined(eta,x);
    end
    toc
    i
end
% map_02 = zeros(101);
% map_11 = zeros(101);
% map_20 = zeros(101);
% map = zeros(101);
% % for i = 1:101
% %     tic
% %     for j = 1:101
% %         map(i,j) = Ns*NLL_S2([pi/100*(i-1);pi/100*(j-1)],x);
% %     end
% %     toc
% %     i
% % end
% % load('testChi.mat');
% % chi = 0;
% for i = 1:101
%     tic
%     for j = 1:101
%         map_02(i,j) = Ns_02*NLL_S3([pi/100*(i-1);pi/100*(j-1);pi/100*(i-1)],x_02);
%     end
%     toc
%     i
% end
% for i = 1:101
%     tic
%     for j = 1:101
%         map_11(i,j) = Ns_11*NLL_S3([pi/100*(i-1);pi/100*(j-1);pi/100*(i-1)],x_11);
%     end
%     toc
%     i
% end
% for i = 1:101
%     tic
%     for j = 1:101
%         map_20(i,j) = Ns_20*NLL_S3([pi/100*(i-1);pi/100*(j-1);pi/100*(i-1)],x_20);
%     end
%     toc
%     i
% end
% for i = 1:101
%     tic
%     for j = 1:101
%         map(i,j) = Ns_02*NLL_S3([pi/100*(i-1);pi/100*(j-1);pi/100*(i-1)],x);
%     end
%     toc
%     i
% end
%
%     0.4194
%     2.4650
%     0.4447
%     0.4188

%     0.5174
%     1.2212
%     0.3947
%     0.5070   
% diary mapS4
% map_phi = zeros(25);
% for i = 1:25
%     tic
%     for j = 1:25
%         map_phi(i,j) = Ns*NLL_S4([0.25+0.6*(i-1)/24; 1.2212; 0.3947; 0.1 + 0.6/24*(j-1)],x);
%     end
%     toc
%     i
% end
% map_phi
% map_lambda = zeros(25);
% for i = 1:25
%     tic
%     for j = 1:25
%         map_lambda(i,j) = Ns*NLL_S4([0.5174; 1 + 1.8*(i-1)/24; 0.3947; 0.1 + 0.6/24*(j-1)],x);
%     end
%     toc
%     i
% end
% map_lambda
% map_v = zeros(25);
% for i = 1:25
%     tic
%     for j = 1:25
%         map_v(i,j) = Ns*NLL_S4([0.5174; 1.2212; 0.15 + 0.6*(i-1)/24; 0.1 + 0.6/24*(j-1)],x);
%     end
%     toc
%     i
% end
% map_v
% diary off
% save s4_map.mat map_phi map_lambda map_v
% save s3_map_1slice.mat map
%%
load('optimizationMaps.mat');
N = 120;
sig_theta = linspace(0.01, pi-0.01, N);
sig_phi = linspace(0.01, pi-0.01, N);
sig_lam = linspace(0.01, pi-0.01, N);

%%
figure('position',[461 535 235 130]);
cmap=hot(100);
map = map_CL;
%%
[M,cf] = contourf(sig_phi,sig_theta,map.','ShowText','off','EdgeColor','none');
view(2);
cf.LevelList = linspace(min(map,[],'all'),-log(1/2/pi),100);
cf.LineStyle = '--';
colormap(cmap);
colorbar('ticks',[min(map,[],'all'),quantile(map,0.75,'all')]);
title('Contralateral', 'Interpreter','latex');
axis equal;
ax = gca;
ax.XTick = [0  pi];
ax.YTick = [0  pi];
ax.XTickLabel = {'0','$\pi/2$','$\pi$'};
ax.YTickLabel = {'0','$\pi/2$','$\pi$'};
ax.XLabel.String = '$\sigma_\Phi$';
ax.YLabel.String = '$\sigma_\Theta$';
%%
figure('position',[461 535 235 130]);
% cmap=flip(hot(100));
map = map_AP;
[M,cf] = contourf(sig_phi,sig_theta,map.','ShowText','off','EdgeColor','none');
view(2);
cf.LevelList = linspace(min(map,[],'all'),quantile(map,0.75,'all'),100);
cf.LineStyle = '--';
colormap(cmap);
colorbar('ticks',[min(map,[],'all'),quantile(map,0.75,'all')]);
title('Antipodal', 'Interpreter','latex');
axis equal;
ax = gca;
ax.XTick = [0  pi];
ax.YTick = [0  pi];
ax.XTickLabel = {'0','$\pi/2$','$\pi$'};
ax.YTickLabel = {'0','$\pi/2$','$\pi$'};
ax.XLabel.String = '$\sigma_\Phi$';
ax.YLabel.String = '$\sigma_\Theta$';

figure('position',[461 535 235 130]);
% cmap=flip(hot(100));
map = map_PLsym;
[M,cf] = contourf(sig_phi,sig_theta,map.','ShowText','off','EdgeColor','none');
view(2);
cf.LevelList = linspace(min(map,[],'all'),quantile(map,0.9,'all'),100);
cf.LineStyle = '--';
colormap(cmap);
colorbar('ticks',[min(map,[],'all'),quantile(map,0.75,'all')]);
title('Parallel', 'Interpreter','latex');
axis equal;
ax = gca;
ax.XTick = [0  pi];
ax.YTick = [0  pi];
ax.XTickLabel = {'0','$\pi/2$','$\pi$'};
ax.YTickLabel = {'0','$\pi/2$','$\pi$'};
ax.XLabel.String = '$\sigma_\Lambda$';
ax.YLabel.String = '$\sigma_\Phi,\sigma_\Theta$';

figure('position',[461 535 235 130]);
% cmap=flip(hot(100));
map = map_DOsym;
[M,cf] = contourf(sig_phi,sig_theta,map.','ShowText','off','EdgeColor','none');
view(2);
cf.LevelList = linspace(min(map,[],'all'),quantile(map,0.9,'all'),100);
cf.LineStyle = '--';
colormap(cmap);
colorbar('ticks',[min(map,[],'all'),quantile(map,0.75,'all')]);
title('Distance-optimal', 'Interpreter','latex');
axis equal;
ax = gca;
ax.XTick = [0  pi];
ax.YTick = [0  pi];
ax.XTickLabel = {'0','$\pi/2$','$\pi$'};
ax.YTickLabel = {'0','$\pi/2$','$\pi$'};
ax.XLabel.String = '$\sigma_\Lambda$';
ax.YLabel.String = '$\sigma_\Phi,\sigma_\Theta$';

figure('position',[461 535 235 130]);
% cmap=flip(hot(100));
map = map_OTsym;
[M,cf] = contourf(sig_phi,sig_theta,map.','ShowText','off','EdgeColor','none');
view(2);
cf.LevelList = linspace(min(map,[],'all'),quantile(map,0.9,'all'),100);
cf.LineStyle = '--';
colormap(cmap);
colorbar('ticks',[min(map,[],'all'),quantile(map,0.9,'all')]);
title('Orthogonal', 'Interpreter','latex');
axis equal;
ax = gca;
ax.XTick = [0  pi];
ax.YTick = [0  pi];
ax.XTickLabel = {'0','$\pi/2$','$\pi$'};
ax.YTickLabel = {'0','$\pi/2$','$\pi$'};
ax.XLabel.String = '$\sigma_\Lambda$';
ax.YLabel.String = '$\sigma_\Phi,\sigma_\Theta$';


%%
load('BSresultsAll.mat');

figure('position',[461 535 235 130]), hold on;
plot(etaopt_CL_1(1,:),etaopt_CL_1(2,:),'r.');
plot(etaopt_CL_2(1,:),etaopt_CL_2(2,:),'g.');
plot(etaopt_CL_3(1,:),etaopt_CL_3(2,:),'b.');

plot(etaopt_CL(1,:),etaopt_CL(2,:),'k.');
axis equal;
xlim([0,pi]);
ylim([0,pi]);
xticks([0,pi]);
yticks([0,pi]);
title('contralateral');

figure('position',[461 535 235 130]), hold on;
plot(etaopt_AP_1(1,:),etaopt_AP_1(2,:),'r.');
plot(etaopt_AP_2(1,:),etaopt_AP_2(2,:),'g.');
plot(etaopt_AP_3(1,:),etaopt_AP_3(2,:),'b.');

plot(etaopt_AP(1,:),etaopt_AP(2,:),'k.');
axis equal;
xlim([0,pi]);
ylim([0,pi]);
xticks([0,pi]);
yticks([0,pi]);
title('antipodal');

figure('position',[461 535 235 130]), hold on;
plot(etaopt_PL_1(2,:),etaopt_PL_1(3,:),'r.');
plot(etaopt_PL_2(2,:),etaopt_PL_2(3,:),'g.');
plot(etaopt_PL_3(2,:),etaopt_PL_3(3,:),'b.');

plot(etaopt_PL(2,:),etaopt_PL(3,:),'k.');
axis equal;
xlim([0,pi]);
ylim([0,pi]);
xticks([0,pi]);
yticks([0,pi]);
title('parallel');


figure('position',[461 535 235 130]), hold on;
plot(etaopt_DO_1(2,:),etaopt_DO_1(3,:),'r.');
plot(etaopt_DO_2(2,:),etaopt_DO_2(3,:),'g.');
plot(etaopt_DO_3(2,:),etaopt_DO_3(3,:),'b.');

plot(etaopt_DO(2,:),etaopt_DO(3,:),'k.');
axis equal;
xlim([0,pi]);
ylim([0,pi]);
xticks([0,pi]);
yticks([0,pi]);
title('distance-optimal');

figure('position',[461 535 235 130]), hold on;
plot(etaopt_OT_1(2,:),etaopt_OT_1(3,:),'r.');
plot(etaopt_OT_2(2,:),etaopt_OT_2(3,:),'g.');
plot(etaopt_OT_3(2,:),etaopt_OT_3(3,:),'b.');

plot(etaopt_OT(2,:),etaopt_OT(3,:),'k.');
axis equal;
xlim([0,pi]);
ylim([0,pi]);
xticks([0,pi]);
yticks([0,pi]);
title('orthogonal');

%% By dataset
load('BSresultsAll.mat');

N = length(nllopt_CL);
mAICset1 = zeros(5,N);
mAICset2 = zeros(5,N);
mAICset3 = zeros(5,N);
mAICsetAll = zeros(5,N);
sAICset1 = zeros(5,N);
sAICset2 = zeros(5,N);
sAICset3 = zeros(5,N);
sAICsetAll = zeros(5,N);

Datasize = [251;233;215;699];
AICopt_CL = nllopt_CL.*Datasize*2 + 2;
AICopt_AP = nllopt_AP.*Datasize*2 + 2;
AICopt_PL = nllopt_PL.*Datasize*2 + 4;
AICopt_DO = nllopt_DO.*Datasize*2 + 4;
AICopt_OT = nllopt_OT.*Datasize*2 + 4;

for i = 1:N
    mAICset1(1,i) = mean(AICopt_CL(1,1:i)); % Contralateral
    mAICset1(2,i) = mean(AICopt_AP(1,1:i)); % Antipodal
    mAICset1(3,i) = mean(AICopt_PL(1,1:i)); % Parallel
    mAICset1(4,i) = mean(AICopt_DO(1,1:i)); % Distance-optimal
    mAICset1(5,i) = mean(AICopt_OT(1,1:i)); % Orthogonal
    
    mAICset2(1,i) = mean(AICopt_CL(2,1:i));
    mAICset2(2,i) = mean(AICopt_AP(2,1:i));
    mAICset2(3,i) = mean(AICopt_PL(2,1:i));
    mAICset2(4,i) = mean(AICopt_DO(2,1:i));
    mAICset2(5,i) = mean(AICopt_OT(2,1:i));
    
    mAICset3(1,i) = mean(AICopt_CL(3,1:i));
    mAICset3(2,i) = mean(AICopt_AP(3,1:i));
    mAICset3(3,i) = mean(AICopt_PL(3,1:i));
    mAICset3(4,i) = mean(AICopt_DO(3,1:i));
    mAICset3(5,i) = mean(AICopt_OT(3,1:i));
    
    mAICsetAll(1,i) = mean(AICopt_CL(4,1:i));
    mAICsetAll(2,i) = mean(AICopt_AP(4,1:i));
    mAICsetAll(3,i) = mean(AICopt_PL(4,1:i));
    mAICsetAll(4,i) = mean(AICopt_DO(4,1:i));
    mAICsetAll(5,i) = mean(AICopt_OT(4,1:i));
    
    sAICset1(1,i) = std(AICopt_CL(1,1:i)); % Contralateral
    sAICset1(2,i) = std(AICopt_AP(1,1:i)); % Antipodal
    sAICset1(3,i) = std(AICopt_PL(1,1:i)); % Parallel
    sAICset1(4,i) = std(AICopt_DO(1,1:i)); % Distance-optimal
    sAICset1(5,i) = std(AICopt_OT(1,1:i)); % Orthogonal
    
    sAICset2(1,i) = std(AICopt_CL(2,1:i));
    sAICset2(2,i) = std(AICopt_AP(2,1:i));
    sAICset2(3,i) = std(AICopt_PL(2,1:i));
    sAICset2(4,i) = std(AICopt_DO(2,1:i));
    sAICset2(5,i) = std(AICopt_OT(2,1:i));
    
    sAICset3(1,i) = std(AICopt_CL(3,1:i));
    sAICset3(2,i) = std(AICopt_AP(3,1:i));
    sAICset3(3,i) = std(AICopt_PL(3,1:i));
    sAICset3(4,i) = std(AICopt_DO(3,1:i));
    sAICset3(5,i) = std(AICopt_OT(3,1:i));
    
    sAICsetAll(1,i) = std(AICopt_CL(4,1:i));
    sAICsetAll(2,i) = std(AICopt_AP(4,1:i));
    sAICsetAll(3,i) = std(AICopt_PL(4,1:i));
    sAICsetAll(4,i) = std(AICopt_DO(4,1:i));
    sAICsetAll(5,i) = std(AICopt_OT(4,1:i));
end
AICm_norm1 = mAICset1(:,end)'/Datasize(1);
AICm_norm2 = mAICset2(:,end)'/Datasize(2);
AICm_norm3 = mAICset3(:,end)'/Datasize(3);
AICm_normAll = mAICsetAll(:,end)'/Datasize(4);
AICs_norm1 = sAICset1(:,end)'/Datasize(1);
AICs_norm2 = sAICset2(:,end)'/Datasize(2);
AICs_norm3 = sAICset3(:,end)'/Datasize(3);
AICs_normAll = sAICsetAll(:,end)'/Datasize(4);
% slow predator
shade_g = figure();shade_g.Position = [440 678 390 120];
subplot(1,4,1),x_axis = 1:length(mAICset1);lColor = [0.2 0.6 0.7];aColor = [0.5 0.75 0.85];lWidth = 1.5;
[patch_CL,l_CL] = plot_shadederrorbar(mAICset1(1,:),sAICset1(1,:),x_axis,lColor,aColor,lWidth);
lColor = [0.65 0.3 0.75];aColor = [0.85 0.6 0.95];
[patch_AP,l_AP] = plot_shadederrorbar(mAICset1(2,:),sAICset1(2,:),x_axis,lColor,aColor,lWidth);
lColor = [0.9 0.4 0.1];aColor = [1 0.7 0.5];
[patch_PL,l_PL] = plot_shadederrorbar(mAICset1(3,:),sAICset1(3,:),x_axis,lColor,aColor,lWidth);
lColor = [0 143 0]./255;aColor = [106 243 66]./255;
[patch_DO,l_DO] = plot_shadederrorbar(mAICset1(4,:),sAICset1(4,:),x_axis,lColor,aColor,lWidth);
lColor = [0 0 0];aColor = [0.6 0.6 0.6];
[patch_OT,l_OT] = plot_shadederrorbar(mAICset1(5,:),sAICset1(5,:),x_axis,lColor,aColor,lWidth);
title('dataset 1');xlabel('bootstrapping');ylabel('AIC');
legend([l_CL,l_AP,l_PL,l_DO,l_OT],...
    'contralateral','antipodal', 'parallel','distance-optimal','orthogonal','Location','bestoutside');
% 
% 
% 
% % dataset 2
% subplot(1,4,2),x_axis = 1:length(mnll_s1);lColor = [0.2 0.6 0.7];aColor = [0.5 0.75 0.85];lWidth = 1.5;
% [patch_s1,l_s1] = plot_shadederrorbar(mnll_s1(2,:),snll_s1(2,:),x_axis,lColor,aColor,lWidth);
% x_axis = 1:length(mnll_s2);lColor = [0.65 0.3 0.75];aColor = [0.85 0.6 0.95];
% [patch_s2,l_s2] = plot_shadederrorbar(mnll_s2(2,:),snll_s2(2,:),x_axis,lColor,aColor,lWidth);
% x_axis = 1:length(mnll_s3);lColor = [0.9 0.4 0.1];aColor = [1 0.7 0.5];
% [patch_s3,l_s3] = plot_shadederrorbar(mnll_s3(2,:),snll_s3(2,:),x_axis,lColor,aColor,lWidth);
% x_axis = 1:length(mnll_s4_1);lColor = [0 143 0]./255;aColor = [106 243 66]./255;
% [patch_s4_1,l_s4_1] = plot_shadederrorbar(mnll_s4_1(2,:),snll_s4_1(2,:),x_axis,lColor,aColor,lWidth);
% x_axis = 1:length(mnll_s4_2);lColor = [0 0 0];aColor = [0.6 0.6 0.6];
% % [patch_s4_2,l_s4_2] = plot_shadederrorbar(mnll_s4_2(2,:),snll_s4_2(2,:),x_axis,lColor,aColor,lWidth);
% % x_axis = 1:length(mnll_s4_3);lColor = [0.7 0.7 0.1];aColor = [1 1 0.2];
% [patch_s4_3,l_s4_3] = plot_shadederrorbar(mnll_s4_3(2,:),snll_s4_3(2,:),x_axis,lColor,aColor,lWidth);
% title('dataset 2');xlabel('bootstrapping');ylabel('NLL');
% % legend([l_s1,l_s2,l_s3,l_s4_1,l_s4_2,l_s4_3],'1','2', '3','4-1','4-2','4-3','Location','bestoutside');
% 
% 
% 
% 
% subplot(1,4,3),x_axis = 1:length(mnll_s1);lColor = [0.2 0.6 0.7];aColor = [0.5 0.75 0.85];lWidth = 1.5;
% [patch_s1,l_s1] = plot_shadederrorbar(mnll_s1(3,:),snll_s1(3,:),x_axis,lColor,aColor,lWidth);
% x_axis = 1:length(mnll_s2);lColor = [0.65 0.3 0.75];aColor = [0.85 0.6 0.95];
% [patch_s2,l_s2] = plot_shadederrorbar(mnll_s2(3,:),snll_s2(3,:),x_axis,lColor,aColor,lWidth);
% x_axis = 1:length(mnll_s3);lColor = [0.9 0.4 0.1];aColor = [1 0.7 0.5];
% [patch_s3,l_s3] = plot_shadederrorbar(mnll_s3(3,:),snll_s3(3,:),x_axis,lColor,aColor,lWidth);
% x_axis = 1:length(mnll_s4_1);lColor = [0 143 0]./255;aColor = [106 243 66]./255;
% [patch_s4_1,l_s4_1] = plot_shadederrorbar(mnll_s4_1(3,:),snll_s4_1(3,:),x_axis,lColor,aColor,lWidth);
% x_axis = 1:length(mnll_s4_2);lColor = [0 0 0];aColor = [0.6 0.6 0.6];
% % [patch_s4_2,l_s4_2] = plot_shadederrorbar(mnll_s4_2(3,:),snll_s4_2(3,:),x_axis,lColor,aColor,lWidth);
% % x_axis = 1:length(mnll_s4_3);lColor = [0.7 0.7 0.1];aColor = [1 1 0.2];
% [patch_s4_3,l_s4_3] = plot_shadederrorbar(mnll_s4_3(3,:),snll_s4_3(3,:),x_axis,lColor,aColor,lWidth);
% title('dataset 3');xlabel('bootstrapping');ylabel('NLL');
% % legend([l_s1,l_s2,l_s3,l_s4_1,l_s4_2,l_s4_3],'1','2', '3','4-1','4-2','4-3','Location','bestoutside');
% 
% 
% subplot(1,4,4),x_axis = 1:length(mnll_s1);lColor = [0.2 0.6 0.7];aColor = [0.5 0.75 0.85];lWidth = 1.5;
% [patch_s1,l_s1] = plot_shadederrorbar(mnll_s1(4,:),snll_s1(4,:),x_axis,lColor,aColor,lWidth);
% x_axis = 1:length(mnll_s2);lColor = [0.65 0.3 0.75];aColor = [0.85 0.6 0.95];
% [patch_s2,l_s2] = plot_shadederrorbar(mnll_s2(4,:),snll_s2(4,:),x_axis,lColor,aColor,lWidth);
% x_axis = 1:length(mnll_s3);lColor = [0.9 0.4 0.1];aColor = [1 0.7 0.5];
% [patch_s3,l_s3] = plot_shadederrorbar(mnll_s3(4,:),snll_s3(4,:),x_axis,lColor,aColor,lWidth);
% x_axis = 1:length(mnll_s4_1);lColor = [0 143 0]./255;aColor = [106 243 66]./255;
% [patch_s4_1,l_s4_1] = plot_shadederrorbar(mnll_s4_1(4,:),snll_s4_1(4,:),x_axis,lColor,aColor,lWidth);
% x_axis = 1:length(mnll_s4_2);lColor = [0 0 0];aColor = [0.6 0.6 0.6];
% % [patch_s4_2,l_s4_2] = plot_shadederrorbar(mnll_s4_2(4,:),snll_s4_2(4,:),x_axis,lColor,aColor,lWidth);
% % x_axis = 1:length(mnll_s4_3);lColor = [0.7 0.7 0.1];aColor = [1 1 0.2];
% [patch_s4_3,l_s4_3] = plot_shadederrorbar(mnll_s4_3(4,:),snll_s4_3(4,:),x_axis,lColor,aColor,lWidth);
% title('dataset 4(all)');xlabel('bootstrapping');ylabel('NLL');
% % legend([l_s1,l_s2,l_s3,l_s4_1,l_s4_2,l_s4_3],'1','2', '3','4-1','4-2','4-3','Location','bestoutside');