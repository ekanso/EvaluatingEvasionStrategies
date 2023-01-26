clear;
close all;

set(groot,'defaultLineLineWidth',2);
set(groot,'defaultFigureColor','w');
set(groot,'defaultTextFontsize',12);
set(groot,'defaultAxesFontsize',12);
set(groot,'defaultPolarAxesFontsize',12);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultPolarAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesLineWidth',1);
colorcode_strategy = [239,65,56;
    41,152,213;
    75,190,108;
    150,109,183;
    255,143,41]/255;
%% new version figure
% load('distance-optimal(0.5)BS.mat');
% load('parallelBS.mat');
% load('BSresultsAll.mat');

N = length(nllopt_CL);
mnllset1 = zeros(5,N);
mnllset2 = zeros(5,N);
mnllset3 = zeros(5,N);
mnllsetAll = zeros(5,N);
snllset1 = zeros(5,N);
snllset2 = zeros(5,N);
snllset3 = zeros(5,N);
snllsetAll = zeros(5,N);

for i = 1:N
    mnllset1(1,i) = mean(nllopt_CL(1,1:i)); % Contralateral
    mnllset1(2,i) = mean(nllopt_AP(1,1:i)); % Antipodal
    mnllset1(3,i) = mean(nllopt_PL(1,1:i)); % Parallel
    mnllset1(4,i) = mean(nllopt_DO(1,1:i)); % Distance-optimal
    mnllset1(5,i) = mean(nllopt_OT(1,1:i)); % Orthogonal
    
    mnllset2(1,i) = mean(nllopt_CL(2,1:i));
    mnllset2(2,i) = mean(nllopt_AP(2,1:i));
    mnllset2(3,i) = mean(nllopt_PL(2,1:i));
    mnllset2(4,i) = mean(nllopt_DO(2,1:i));
    mnllset2(5,i) = mean(nllopt_OT(2,1:i));
    
    mnllset3(1,i) = mean(nllopt_CL(3,1:i));
    mnllset3(2,i) = mean(nllopt_AP(3,1:i));
    mnllset3(3,i) = mean(nllopt_PL(3,1:i));
    mnllset3(4,i) = mean(nllopt_DO(3,1:i));
    mnllset3(5,i) = mean(nllopt_OT(3,1:i));
    
    mnllsetAll(1,i) = mean(nllopt_CL(4,1:i));
    mnllsetAll(2,i) = mean(nllopt_AP(4,1:i));
    mnllsetAll(3,i) = mean(nllopt_PL(4,1:i));
    mnllsetAll(4,i) = mean(nllopt_DO(4,1:i));
    mnllsetAll(5,i) = mean(nllopt_OT(4,1:i));
    
    snllset1(1,i) = std(nllopt_CL(1,1:i)); % Contralateral
    snllset1(2,i) = std(nllopt_AP(1,1:i)); % Antipodal
    snllset1(3,i) = std(nllopt_PL(1,1:i)); % Parallel
    snllset1(4,i) = std(nllopt_DO(1,1:i)); % Distance-optimal
    snllset1(5,i) = std(nllopt_OT(1,1:i)); % Orthogonal
    
    snllset2(1,i) = std(nllopt_CL(2,1:i));
    snllset2(2,i) = std(nllopt_AP(2,1:i));
    snllset2(3,i) = std(nllopt_PL(2,1:i));
    snllset2(4,i) = std(nllopt_DO(2,1:i));
    snllset2(5,i) = std(nllopt_OT(2,1:i));
    
    snllset3(1,i) = std(nllopt_CL(3,1:i));
    snllset3(2,i) = std(nllopt_AP(3,1:i));
    snllset3(3,i) = std(nllopt_PL(3,1:i));
    snllset3(4,i) = std(nllopt_DO(3,1:i));
    snllset3(5,i) = std(nllopt_OT(3,1:i));
    
    snllsetAll(1,i) = std(nllopt_CL(4,1:i));
    snllsetAll(2,i) = std(nllopt_AP(4,1:i));
    snllsetAll(3,i) = std(nllopt_PL(4,1:i));
    snllsetAll(4,i) = std(nllopt_DO(4,1:i));
    snllsetAll(5,i) = std(nllopt_OT(4,1:i));
end


x = [1:5 , 9:13, 17:21, 25:29];
meanArray = [mnllset1(:,end)', mnllset2(:,end)', mnllset3(:,end)', mnllsetAll(:,end)'];
stdArray = [snllset1(:,end)', snllset2(:,end)', snllset3(:,end)', snllsetAll(:,end)'];

figure('position',[253 1177 321 182]);
eb = errorbar(x, meanArray,stdArray,'marker','s','linestyle','none','markersize',4,'capsize',0);
box off;
% xticks(2.5:8:30); 
xlabel('dataset');
ylim([1.35,1.85]); ylabel('Normalized NLL');
yticks(1.35:0.1:1.85);
%% delta_AIC

% load('BSresultsAll.mat');

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
AICopt_CL = nllopt_CL.*Datasize*2 + 2*2;
AICopt_AP = nllopt_AP.*Datasize*2 + 2*2;
AICopt_PL = nllopt_PL.*Datasize*2 + 3*2;
AICopt_DO = nllopt_DO.*Datasize*2 + 3*2;
AICopt_OT = nllopt_OT.*Datasize*2 + 3*2;

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


x = [1:5 , 9:13, 17:21, 25:29];
meanArray = [AICm_norm1-min(AICm_norm1), AICm_norm2-min(AICm_norm2),...
    AICm_norm3-min(AICm_norm3), AICm_normAll-min(AICm_normAll)];
stdArray = [AICs_norm1, AICs_norm2, AICs_norm3, AICs_normAll];

figure('position',[253 1177 321 182]);
eb = errorbar(x, meanArray,stdArray,'marker','s','linestyle','none','markersize',4,'capsize',0);
box off;
xticks(3:8:30); 
xlabel('dataset');
ylim([-0.2,0.8]);
ylabel('Normalized $\Delta$AIC');
yticks(-0.2:0.2:0.8);
%% rankings

[lowest_slow,indexes_slow] = min([nllopt_DO(1,:);
    nllopt_OT(1,:);
    nllopt_PL(1,:);
    nllopt_AP(1,:);
    nllopt_CL(1,:);]);
[lowest_mid,indexes_mid] = min([nllopt_DO(2,:);
    nllopt_OT(2,:);
    nllopt_PL(2,:);
    nllopt_AP(2,:);
    nllopt_CL(2,:);]);
[lowest_fast,indexes_fast] = min([nllopt_DO(3,:);
    nllopt_OT(3,:);
    nllopt_PL(3,:);
    nllopt_AP(3,:);
    nllopt_CL(3,:);]);
[lowest_combined,indexes_combined] = min([nllopt_DO(4,:);
    nllopt_OT(4,:);
    nllopt_PL(4,:);
    nllopt_AP(4,:);
    nllopt_CL(4,:);]);
histcounts(indexes_slow,0.5:5.5)
histcounts(indexes_mid,0.5:5.5)
histcounts(indexes_fast,0.5:5.5)
histcounts(indexes_combined,0.5:5.5)

%% exponential plot
% V=2cm/s
AIC_mean = mAICset1(:,end)';
AIC_delta = AIC_mean - min(AIC_mean);
AIC_std = sAICset1(:,end)';
y = -0.1:-0.1:-0.5;
figure('Position',[440 684 348 113]);
eb = errorbar(AIC_delta,y,AIC_std,'horizontal','marker','s','linestyle','none','markersize',4,'capsize',0);
hold on;
x = 0:0.1:200;
plot(x, exp(-x/2));
plot([-50,200],[0,0],'k--');
xlim([-25,160]);
ylim([-0.6,1.2]);
xlabel('$\Delta$AIC');
ylabel('$\exp(-\Delta\rm{AIC})$');
title('slow')
% V=11cm/s
AIC_mean = mAICset2(:,end)';
AIC_delta = AIC_mean - min(AIC_mean);
AIC_std = sAICset2(:,end)';
y = -0.1:-0.1:-0.5;
figure('Position',[440 684 348 113]);
eb = errorbar(AIC_delta,y,AIC_std,'horizontal','marker','s','linestyle','none','markersize',4,'capsize',0);
hold on;
x = 0:0.1:200;
plot(x, exp(-x/2));
plot([-50,200],[0,0],'k--');
xlim([-25,160]);
ylim([-0.6,1.2]);
xlabel('$\Delta$AIC');
ylabel('$\exp(-\Delta\rm{AIC})$');
title('mid-speed')
% V=20cm/s
AIC_mean = mAICset3(:,end)';
AIC_delta = AIC_mean - min(AIC_mean);
AIC_std = sAICset3(:,end)';
y = -0.1:-0.1:-0.5;
figure('Position',[440 684 348 113]);
eb = errorbar(AIC_delta,y,AIC_std,'horizontal','marker','s','linestyle','none','markersize',4,'capsize',0);
hold on;
x = 0:0.1:200;
plot(x, exp(-x/2));
plot([-50,200],[0,0],'k--');
xlim([-25,160]);
ylim([-0.6,1.2]);
xlabel('$\Delta$AIC');
ylabel('$\exp(-\Delta\rm{AIC})$');
title('fast')
% V=2,11,20cm/s
AIC_mean = mAICsetAll(:,end)';
AIC_delta = AIC_mean - min(AIC_mean);
AIC_std = sAICsetAll(:,end)';
y = -0.1:-0.1:-0.5;
figure('Position',[440 684 348 113]);
eb = errorbar(AIC_delta,y,AIC_std,'horizontal','marker','s','linestyle','none','markersize',4,'capsize',0);
hold on;
x = 0:0.1:500;
plot(x, exp(-x/2));
plot([-50,500],[0,0],'k--');
xlim([-40,450]);
ylim([-0.6,1.2]);
xlabel('$\Delta$AIC');
ylabel('$\exp(-\Delta\rm{AIC})$');
title('combined')
%% Convergence, by strategy


fig_convergence = figure('Position',[200 150 2200 1100]);
subplot(2,3,1),x_axis = 1:length(mnllset1);lColor = [52 148 186]./255;aColor = [128 193 219]./255;lWidth = 2;
[patch_1,l_1] = plot_shadederrorbar(mnllset1(1,:),snllset1(1,:),x_axis,lColor,aColor,lWidth);
x_axis = 1:length(mnllset2);lColor = [117 52 137]./255;aColor = [213 155 240]./255;
[patch_2,l_2] = plot_shadederrorbar(mnllset2(1,:),snllset2(1,:),x_axis,lColor,aColor,lWidth);
x_axis = 1:length(mnllset3);lColor = [201 92 46]./255;aColor = [255 150 92]./255;
[patch_3,l_3] = plot_shadederrorbar(mnllset3(1,:),snllset3(1,:),x_axis,lColor,aColor,lWidth);
x_axis = 1:length(mnllsetAll);lColor = [0 143 0]./255;aColor = [106 243 66]./255;
[patch_all,l_all] = plot_shadederrorbar(mnllsetAll(1,:),snllsetAll(1,:),x_axis,lColor,aColor,lWidth);
ylim([1.3,2]);title('Contralateral Strategy');
lg=legend([l_1,l_2,l_3,l_all],'1','2', '3','all'); lg.Position = [patch_1.Parent.Position(1)+patch_1.Parent.Position(3),...
    patch_1.Parent.Position(2)+patch_1.Parent.Position(4)-lg.Position(4),lg.Position(3),lg.Position(4)];

% f1 = figure();f1.Position = [95 237 1288 420];
% subplot(1,2,1),plot(mnllset1*251,'LineWidth',1.5);
% hold on
% plot(mnllset2*233,'LineWidth',1.5);
% plot(mnllset3*215,'LineWidth',1.5);
% plot(mnllsetAll*699,'LineWidth',1.5);
% legend('1','2','3','all');
% title('Strategy 1,mean');
% subplot(1,2,2),plot(snllset1*251,'LineWidth',1.5);
% hold on
% plot(snllset2*233,'LineWidth',1.5);
% plot(snllset3*215,'LineWidth',1.5);
% plot(snllsetAll*699,'LineWidth',1.5);
% legend('1','2','3','all');
% title('Strategy 1,variance');


subplot(2,3,2),x_axis = 1:length(mnllset1);lColor = [52 148 186]./255;aColor = [128 193 219]./255;lWidth = 2;
[patch_1,l_1] = plot_shadederrorbar(mnllset1(2,:),snllset1(2,:),x_axis,lColor,aColor,lWidth);
x_axis = 1:length(mnllset2);lColor = [117 52 137]./255;aColor = [213 155 240]./255;
[patch_2,l_2] = plot_shadederrorbar(mnllset2(2,:),snllset2(2,:),x_axis,lColor,aColor,lWidth);
x_axis = 1:length(mnllset3);lColor = [201 92 46]./255;aColor = [255 150 92]./255;
[patch_3,l_3] = plot_shadederrorbar(mnllset3(2,:),snllset3(2,:),x_axis,lColor,aColor,lWidth);
x_axis = 1:length(mnllsetAll);lColor = [0 143 0]./255;aColor = [106 243 66]./255;
[patch_all,l_all] = plot_shadederrorbar(mnllsetAll(2,:),snllsetAll(2,:),x_axis,lColor,aColor,lWidth);
ylim([1.3,2]);title('Antipodal Strategy');
lg=legend([l_1,l_2,l_3,l_all],'1','2', '3','all'); lg.Position = [patch_1.Parent.Position(1)+patch_1.Parent.Position(3),...
    patch_1.Parent.Position(2)+patch_1.Parent.Position(4)-lg.Position(4),lg.Position(3),lg.Position(4)];

% f2 = figure();f2.Position = [95 237 1288 420];
% subplot(1,2,1),plot(mnllset1*251,'LineWidth',1.5);
% hold on
% plot(mnllset2*233,'LineWidth',1.5);
% plot(mnllset3*215,'LineWidth',1.5);
% plot(mnllsetAll*699,'LineWidth',1.5);
% legend('1','2','3','all');
% title('Strategy 2,mean');
% subplot(1,2,2),plot(snllset1*251,'LineWidth',1.5);
% hold on
% plot(snllset2*233,'LineWidth',1.5);
% plot(snllset3*215,'LineWidth',1.5);
% plot(snllsetAll*699,'LineWidth',1.5);
% legend('1','2','3','all');
% title('Strategy 2,variance');

subplot(2,3,3),x_axis = 1:length(mnllset1);lColor = [52 148 186]./255;aColor = [128 193 219]./255;lWidth = 2;
[patch_1,l_1] = plot_shadederrorbar(mnllset1(3,:),snllset1(3,:),x_axis,lColor,aColor,lWidth);
x_axis = 1:length(mnllset2);lColor = [117 52 137]./255;aColor = [213 155 240]./255;
[patch_2,l_2] = plot_shadederrorbar(mnllset2(3,:),snllset2(3,:),x_axis,lColor,aColor,lWidth);
x_axis = 1:length(mnllset3);lColor = [201 92 46]./255;aColor = [255 150 92]./255;
[patch_3,l_3] = plot_shadederrorbar(mnllset3(3,:),snllset3(3,:),x_axis,lColor,aColor,lWidth);
x_axis = 1:length(mnllsetAll);lColor = [0 143 0]./255;aColor = [106 243 66]./255;
[patch_all,l_all] = plot_shadederrorbar(mnllsetAll(3,:),snllsetAll(3,:),x_axis,lColor,aColor,lWidth);
xlim([0,length(mnllset1)]);ylim([1.3,2]);title('Parallel Strategy');
lg=legend([l_1,l_2,l_3,l_all],'1','2', '3','all'); lg.Position = [patch_1.Parent.Position(1)+patch_1.Parent.Position(3),...
    patch_1.Parent.Position(2)+patch_1.Parent.Position(4)-lg.Position(4),lg.Position(3),lg.Position(4)];

% f3 = figure();f3.Position = [95 237 1288 420];
% subplot(1,2,1),plot(mnllset1,'LineWidth',1.5);
% hold on
% plot(mnllset2,'LineWidth',1.5);
% plot(mnllset3,'LineWidth',1.5);
% plot(mnllsetAll,'LineWidth',1.5);
% legend('1','2','3','all');
% title('Strategy 3,mean');
% subplot(1,2,2),plot(snllset1,'LineWidth',1.5);
% hold on
% plot(snllset2,'LineWidth',1.5);
% plot(snllset3,'LineWidth',1.5);
% plot(snllsetAll,'LineWidth',1.5);
% legend('1','2','3','all');
% title('Strategy 3,variance');


subplot(2,3,4),x_axis = 1:length(mnllset1);lColor = [52 148 186]./255;aColor = [128 193 219]./255;lWidth = 2;
[patch_1,l_1] = plot_shadederrorbar(mnllset1(4,:),snllset1(4,:),x_axis,lColor,aColor,lWidth);
x_axis = 1:length(mnllset2);lColor = [117 52 137]./255;aColor = [213 155 240]./255;
[patch_2,l_2] = plot_shadederrorbar(mnllset2(4,:),snllset2(4,:),x_axis,lColor,aColor,lWidth);
x_axis = 1:length(mnllset3);lColor = [201 92 46]./255;aColor = [255 150 92]./255;
[patch_3,l_3] = plot_shadederrorbar(mnllset3(4,:),snllset3(4,:),x_axis,lColor,aColor,lWidth);
x_axis = 1:length(mnllsetAll);lColor = [0 143 0]./255;aColor = [106 243 66]./255;
[patch_all,l_all] = plot_shadederrorbar(mnllsetAll(4,:),snllsetAll(4,:),x_axis,lColor,aColor,lWidth);
xlim([0,length(mnllset1)]);ylim([1.3,2]);title('Distance-optimal Strategy');
lg=legend([l_1,l_2,l_3,l_all],'1','2', '3','all'); lg.Position = [patch_1.Parent.Position(1)+patch_1.Parent.Position(3),...
    patch_1.Parent.Position(2)+patch_1.Parent.Position(4)-lg.Position(4),lg.Position(3),lg.Position(4)];

% f4_1 = figure();f4_1.Position = [95 237 1288 420];
% subplot(1,2,1),plot(mnllset1,'LineWidth',1.5);
% hold on
% plot(mnllset2,'LineWidth',1.5);
% plot(mnllset3,'LineWidth',1.5);
% plot(mnllsetAll,'LineWidth',1.5);
% legend('1','2','3','all');
% title('Strategy 4($\chi = 0$),mean');
% subplot(1,2,2),plot(snllset1,'LineWidth',1.5);
% hold on
% plot(snllset2,'LineWidth',1.5);
% plot(snllset3,'LineWidth',1.5);
% plot(snllsetAll,'LineWidth',1.5);
% legend('1','2','3','all');
% title('Strategy 4($\chi = 0$,variance)');


subplot(2,3,5),x_axis = 1:length(mnllset1);lColor = [52 148 186]./255;aColor = [128 193 219]./255;lWidth = 2;
[patch_1,l_1] = plot_shadederrorbar(mnllset1(5,:),snllset1(5,:),x_axis,lColor,aColor,lWidth);
x_axis = 1:length(mnllset2);lColor = [117 52 137]./255;aColor = [213 155 240]./255;
[patch_2,l_2] = plot_shadederrorbar(mnllset2(5,:),snllset2(5,:),x_axis,lColor,aColor,lWidth);
x_axis = 1:length(mnllset3);lColor = [201 92 46]./255;aColor = [255 150 92]./255;
[patch_3,l_3] = plot_shadederrorbar(mnllset3(5,:),snllset3(5,:),x_axis,lColor,aColor,lWidth);
x_axis = 1:length(mnllsetAll);lColor = [0 143 0]./255;aColor = [106 243 66]./255;
[patch_all,l_all] = plot_shadederrorbar(mnllsetAll(5,:),snllsetAll(5,:),x_axis,lColor,aColor,lWidth);
xlim([0,length(mnllset1)]);ylim([1.3,2]);title('Orthogonal Strategy($\chi = \arccos(1/2)$)');
lg=legend([l_1,l_2,l_3,l_all],'1','2', '3','all'); lg.Position = [patch_1.Parent.Position(1)+patch_1.Parent.Position(3),...
    patch_1.Parent.Position(2)+patch_1.Parent.Position(4)-lg.Position(4),lg.Position(3),lg.Position(4)];

% f4_2 = figure();f4_2.Position = [95 237 1288 420];
% subplot(1,2,1),plot(mnllset1,'LineWidth',1.5);
% hold on
% plot(mnllset2,'LineWidth',1.5);
% plot(mnllset3,'LineWidth',1.5);
% plot(mnllsetAll,'LineWidth',1.5);
% legend('1','2','3','all');
% title('Strategy 4($\chi = \arccos(2/11)$,mean');
% subplot(1,2,2),plot(snllset1,'LineWidth',1.5);
% hold on
% plot(snllset2,'LineWidth',1.5);
% plot(snllset3,'LineWidth',1.5);
% plot(snllsetAll,'LineWidth',1.5);
% legend('1','2','3','all');
% title('Strategy 4($\chi = \arccos(2/11)$,variance)');
%% plot distance-optimal NLLs with varying chi
load('bestchiresults_all.mat');
N = size(nll_bootstrap,2);
mean_set1 = zeros(1,N);
mean_set2 = zeros(1,N);
mean_set3 = zeros(1,N);
mean_setAll = zeros(1,N);
sd_set1 = zeros(1,N);
sd_set2 = zeros(1,N);
sd_set3 = zeros(1,N);
sd_setAll = zeros(1,N);
for i = 1:N    
    mean_setAll(i) = mean(nll_bootstrap(:,i)); 
    sd_setAll(i) = std(nll_bootstrap(:,i)); 
end
load('bestchiresults_02.mat');
for i = 1:N    
    mean_set1(i) = mean(nll_bootstrap(:,i)); 
    sd_set1(i) = std(nll_bootstrap(:,i)); 
end
load('bestchiresults_11.mat');
for i = 1:N    
    mean_set2(i) = mean(nll_bootstrap(:,i)); 
    sd_set2(i) = std(nll_bootstrap(:,i)); % 
end
load('bestchiresults_20.mat');
for i = 1:N    
    mean_set3(i) = mean(nll_bootstrap(:,i)); %
    sd_set3(i) = std(nll_bootstrap(:,i)); %
end
fig_varyChi = figure('Position',[200 150 2200 1100]);

% subplot(4,1,1),
x_axis = (1:size(nll_bootstrap,2))/N*pi/2;lColor = colorcode_strategy(1,:);aColor = 1-(1-colorcode_strategy(1,:))/2;lWidth = 2;
[patch_1,l_1] = plot_shadederrorbar(mean_set1,sd_set1,x_axis,lColor,aColor,lWidth);
% xlim([0,length(mean_set1)]);ylim([1.3,2]);title('slow predator');
% subplot(4,1,2),
x_axis = (1:size(nll_bootstrap,2))/N*pi/2;lColor = colorcode_strategy(2,:);aColor = 1-(1-colorcode_strategy(2,:))/2;lWidth = 2;
[patch_2,l_2] = plot_shadederrorbar(mean_set2,sd_set2,x_axis,lColor,aColor,lWidth);
% xlim([0,length(mean_set1)]);ylim([1.3,2]);title('mid-speed predator');
% subplot(4,1,3),
x_axis = (1:size(nll_bootstrap,2))/N*pi/2;lColor = colorcode_strategy(3,:);aColor = 1-(1-colorcode_strategy(3,:))/2;lWidth = 2;
[patch_3,l_3] = plot_shadederrorbar(mean_set3,sd_set3,x_axis,lColor,aColor,lWidth);
% xlim([0,length(mean_set1)]);ylim([1.3,2]);title('fast predator');
% subplot(4,1,4),
x_axis = (1:size(nll_bootstrap,2))/N*pi/2;lColor = colorcode_strategy(4,:);aColor = 1-(1-colorcode_strategy(4,:))/2;lWidth = 2;
[patch_all,l_all] = plot_shadederrorbar(mean_setAll,sd_setAll,x_axis,lColor,aColor,lWidth);
% xlim([0,length(mean_set1)]);ylim([1.3,2]);title('All data');

% lg=legend([l_1,l_2,l_3,l_all],'1','2', '3','all'); lg.Position = [patch_1.Parent.Position(1)+patch_1.Parent.Position(3),...
%     patch_1.Parent.Position(2)+patch_1.Parent.Position(4)-lg.Position(4),lg.Position(3),lg.Position(4)];
l_1.Marker = 'd';l_2.Marker = '*';l_3.Marker = 'o';l_all.Marker = 's';
l_1.LineStyle = '--';l_2.LineStyle = '--';l_3.LineStyle = '--';
xticks([0,pi/4,pi/2]);
%%
fe1 = figure();fe1.Position = [285,116,841,652];
subplot(2,2,1),scatter(etaopt_s1(1,:),etaopt_s1(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S1,group:all');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,2),scatter(etaopt_s1_1(1,:),etaopt_s1_1(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S1,group:1');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,3),scatter(etaopt_s1_2(1,:),etaopt_s1_2(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S1,group:2');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,4),scatter(etaopt_s1_3(1,:),etaopt_s1_3(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S1,group:3');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');

fe2 = figure();fe2.Position = [285,116,841,652];
subplot(2,2,1),scatter(etaopt_s2(1,:),etaopt_s2(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S2,group:all');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,2),scatter(etaopt_s2_1(1,:),etaopt_s2_1(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S2,group:1');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,3),scatter(etaopt_s2_2(1,:),etaopt_s2_2(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S2,group:2');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,4),scatter(etaopt_s2_3(1,:),etaopt_s2_3(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S2,group:3');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');

fe3pl = figure();fe3pl.Position = [285,116,841,652];
subplot(2,2,1),scatter(etaopt_s3(1,:),etaopt_s3(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S3,group:all');xlabel('$\sigma_\phi$');ylabel('$\sigma_\lambda$');
subplot(2,2,2),scatter(etaopt_s3_1(1,:),etaopt_s3_1(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S3,group:1');xlabel('$\sigma_\phi$');ylabel('$\sigma_\lambda$');
subplot(2,2,3),scatter(etaopt_s3_2(1,:),etaopt_s3_2(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S3,group:2');xlabel('$\sigma_\phi$');ylabel('$\sigma_\lambda$');
subplot(2,2,4),scatter(etaopt_s3_3(1,:),etaopt_s3_3(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S3,group:3');xlabel('$\sigma_\phi$');ylabel('$\sigma_\lambda$');

fe3pt = figure();fe3pt.Position = [285,116,841,652];
subplot(2,2,1),scatter(etaopt_s3(1,:),etaopt_s3(3,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S3,group:all');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,2),scatter(etaopt_s3_1(1,:),etaopt_s3_1(3,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S3,group:1');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,3),scatter(etaopt_s3_2(1,:),etaopt_s3_2(3,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S3,group:2');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,4),scatter(etaopt_s3_3(1,:),etaopt_s3_3(3,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S3,group:3');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');

fe4_1pl = figure();fe4_1pl.Position = [285,116,841,652];
subplot(2,2,1),scatter(etaopt_s4_1(1,:),etaopt_s4_1(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=0$,group:all');xlabel('$\sigma_\phi$');ylabel('$\sigma_\lambda$');
subplot(2,2,2),scatter(etaopt_s4_1_1(1,:),etaopt_s4_1_1(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=0$,group:1');xlabel('$\sigma_\phi$');ylabel('$\sigma_\lambda$');
subplot(2,2,3),scatter(etaopt_s4_1_2(1,:),etaopt_s4_1_2(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=0$,group:2');xlabel('$\sigma_\phi$');ylabel('$\sigma_\lambda$');
subplot(2,2,4),scatter(etaopt_s4_1_3(1,:),etaopt_s4_1_3(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=0$,group:3');xlabel('$\sigma_\phi$');ylabel('$\sigma_\lambda$');

fe4_1pt = figure();fe4_1pt.Position = [285,116,841,652];
subplot(2,2,1),scatter(etaopt_s4_1(1,:),etaopt_s4_1(3,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=0$,group:all');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,2),scatter(etaopt_s4_1_1(1,:),etaopt_s4_1_1(3,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=0$,group:1');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,3),scatter(etaopt_s4_1_2(1,:),etaopt_s4_1_2(3,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=0$,group:2');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,4),scatter(etaopt_s4_1_3(1,:),etaopt_s4_1_3(3,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=0$,group:3');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');

fe4_2pl = figure();fe4_2pl.Position = [285,116,841,652];
subplot(2,2,1),scatter(etaopt_s4_2(1,:),etaopt_s4_2(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=\arccos(2/11)$,group:all');xlabel('$\sigma_\phi$');ylabel('$\sigma_\lambda$');
subplot(2,2,2),scatter(etaopt_s4_2_1(1,:),etaopt_s4_2_1(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=\arccos(2/11)$,group:1');xlabel('$\sigma_\phi$');ylabel('$\sigma_\lambda$');
subplot(2,2,3),scatter(etaopt_s4_2_2(1,:),etaopt_s4_2_2(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=\arccos(2/11)$,group:2');xlabel('$\sigma_\phi$');ylabel('$\sigma_\lambda$');
subplot(2,2,4),scatter(etaopt_s4_2_3(1,:),etaopt_s4_2_3(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=\arccos(2/11)$,group:3');xlabel('$\sigma_\phi$');ylabel('$\sigma_\lambda$');

fe4_2pt = figure();fe4_2pt.Position = [285,116,841,652];
subplot(2,2,1),scatter(etaopt_s4_2(1,:),etaopt_s4_2(3,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=\arccos(2/11)$,group:all');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,2),scatter(etaopt_s4_2_1(1,:),etaopt_s4_2_1(3,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=\arccos(2/11)$,group:1');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,3),scatter(etaopt_s4_2_2(1,:),etaopt_s4_2_2(3,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=\arccos(2/11)$,group:2');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,4),scatter(etaopt_s4_2_3(1,:),etaopt_s4_2_3(3,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=\arccos(2/11)$,group:3');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');

fe4_3pl = figure();fe4_3pl.Position = [285,116,841,652];
subplot(2,2,1),scatter(etaopt_s4_3(1,:),etaopt_s4_3(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=\arccos(2/20)$,group:all');xlabel('$\sigma_\phi$');ylabel('$\sigma_\lambda$');
subplot(2,2,2),scatter(etaopt_s4_3_1(1,:),etaopt_s4_3_1(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=\arccos(2/20)$,group:1');xlabel('$\sigma_\phi$');ylabel('$\sigma_\lambda$');
subplot(2,2,3),scatter(etaopt_s4_3_2(1,:),etaopt_s4_3_2(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=\arccos(2/20)$,group:2');xlabel('$\sigma_\phi$');ylabel('$\sigma_\lambda$');
subplot(2,2,4),scatter(etaopt_s4_3_3(1,:),etaopt_s4_3_3(2,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=\arccos(2/20)$,group:3');xlabel('$\sigma_\phi$');ylabel('$\sigma_\lambda$');

fe4_3pt = figure();fe4_3pt.Position = [285,116,841,652];
subplot(2,2,1),scatter(etaopt_s4_3(1,:),etaopt_s4_3(3,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=\arccos(2/20)$,group:all');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,2),scatter(etaopt_s4_3_1(1,:),etaopt_s4_3_1(3,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=\arccos(2/20)$,group:1');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,3),scatter(etaopt_s4_3_2(1,:),etaopt_s4_3_2(3,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=\arccos(2/20)$,group:2');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');
subplot(2,2,4),scatter(etaopt_s4_3_3(1,:),etaopt_s4_3_3(3,:),72,'.');
xlim([0,pi/2]);ylim([0,pi/2]);axis equal;
title('S4,$\chi=\arccos(2/20)$,group:3');xlabel('$\sigma_\phi$');ylabel('$\sigma_\theta$');


%%
figure,subplot(2,3,1),histogram(251*nllopt_CL(1,:),10),title('Strategy 1');
subplot(2,3,2),histogram(251*nllopt_AP(1,:),10),title('Strategy 2');
subplot(2,3,3),histogram(251*nllopt_OT(1,:),10),title('Strategy 3');
subplot(2,3,4),histogram(251*nllopt_PL(1,:),10),title('Strategy 4-1');
subplot(2,3,5),histogram(251*nllopt_DO(1,:),10),title('Strategy 4-2');
subplot(2,3,6),histogram(251*nllopt_s4_3(1,:),10),title('Strategy 4-3');

figure,subplot(2,3,1),histogram(233*nllopt_CL(2,:),10),title('Strategy 1');
subplot(2,3,2),histogram(233*nllopt_AP(2,:),10),title('Strategy 2');
subplot(2,3,3),histogram(233*nllopt_OT(2,:),10),title('Strategy 3');
subplot(2,3,4),histogram(233*nllopt_PL(2,:),10),title('Strategy 4-1');
subplot(2,3,5),histogram(233*nllopt_DO(2,:),10),title('Strategy 4-2');
subplot(2,3,6),histogram(233*nllopt_s4_3(2,:),10),title('Strategy 4-3');

figure,subplot(2,3,1),histogram(251*nllopt_CL(3,:),10),title('Strategy 1');
subplot(2,3,2),histogram(215*nllopt_AP(3,:),10),title('Strategy 2');
subplot(2,3,3),histogram(215*nllopt_OT(3,:),10),title('Strategy 3');
subplot(2,3,4),histogram(215*nllopt_PL(3,:),10),title('Strategy 4-1');
subplot(2,3,5),histogram(215*nllopt_DO(3,:),10),title('Strategy 4-2');
subplot(2,3,6),histogram(215*nllopt_s4_3(3,:),10),title('Strategy 4-3');

figure,subplot(2,3,1),histogram(699*nllopt_CL(4,:),10),title('Strategy 1');
subplot(2,3,2),histogram(699*nllopt_AP(4,:),10),title('Strategy 2');
subplot(2,3,3),histogram(699*nllopt_OT(4,:),10),title('Strategy 3');
subplot(2,3,4),histogram(699*nllopt_PL(4,:),10),title('Strategy 4-1');
subplot(2,3,5),histogram(699*nllopt_DO(4,:),10),title('Strategy 4-2');
subplot(2,3,6),histogram(699*nllopt_s4_3(4,:),10),title('Strategy 4-3');


%% error bars
load('optNLLs.mat');
f_e = figure();
nll02 = [nll02(:,1:4) nll02(:,end)];
nll11 = [nll11(:,1:4) nll11(:,end)];
nll20 = [nll20(:,1:4) nll20(:,end)];
nllall = [nllall(:,1:4) nllall(:,end)];
subplot(1,4,1),e02 = errorbar(nll02(1,:),nll02(2,:),'s');
xlim([0.5,5.5]);ylim([360,460]);xticks(1:5);yticks(360:20:460);box off;
subplot(1,4,2),e11 = errorbar(nll11(1,:),nll11(2,:),'s');
xlim([0.5,5.5]);ylim([330,430]);xticks(1:5);yticks(330:20:430);box off;
subplot(1,4,3),e20 = errorbar(nll20(1,:),nll20(2,:),'s');
xlim([0.5,5.5]);ylim([300,400]);xticks(1:5);yticks(300:20:400);box off;
subplot(1,4,4),eall = errorbar(nllall(1,:),nllall(2,:),'s');
xlim([0.5,5.5]);ylim([1000 1300]);xticks(1:5);yticks(1000:60:1300);box off;
e02.MarkerSize = 10;e02.LineWidth = 1.5;
e11.MarkerSize = 10;e11.LineWidth = 1.5;
e20.MarkerSize = 10;e20.LineWidth = 1.5;
eall.MarkerSize = 10;eall.LineWidth = 1.5;