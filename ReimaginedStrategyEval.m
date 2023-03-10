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
load('evasionData.mat');
% u is pre-computed larvae average speed using u = \delta d/\delta t
% v is predator speed
% phi is predator angular position
% lambda is 'relative heading' of predator
% theta is the evasion direction
phi = data.Phi(isfinite(data.U));
u   = data.U(isfinite(data.U));
v   = data.V(isfinite(data.U))*0.01; % original unit is cm/s, converted to m/s
lam = data.Lambda(isfinite(data.U));
theta = data.Delta(isfinite(data.U));

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
theta_PL = phi + lam + pi;
theta_PL= wrapToPi(theta_PL);
%% Distance-optimal, u/v = 0.5
branch = sign(lam);
chi = acos(0.5);
theta_DO = phi + lam + pi - branch*chi;

theta_DO= wrapToPi(theta_DO);
%% Orthogonal
branch = sign(lam);
chi = acos(0);
theta_OT = phi + lam + pi - branch*chi;

theta_OT= wrapToPi(theta_OT);
%% Distance-optimal, u/v from data (real)
branch = sign(lam);
chi = real(acos(u./v));
theta_DOR = phi + lam + pi - branch.*chi;
theta_DOR= wrapTo2Pi(theta_DOR);
%% Re-imagination of strategies
% power
% theta_CL = theta_CL.*(5/9).^(abs(theta_CL)/pi);
% theta_AP = theta_AP.*(5/9).^(abs(theta_AP)/pi);
% theta_PL = theta_PL.*(5/9).^(abs(theta_PL)/pi);
% theta_DO = theta_DO.*(5/9).^(abs(theta_DO)/pi);
% theta_OT = theta_OT.*(5/9).^(abs(theta_OT)/pi);
% polynomial (currently chosen)
theta_CL = theta_CL.*(1 - abs(theta_CL)/pi*4/9);
theta_AP = theta_AP.*(1 - abs(theta_AP)/pi*4/9);
theta_PL = theta_PL.*(1 - abs(theta_PL)/pi*4/9);
theta_DO = theta_DO.*(1 - abs(theta_DO)/pi*4/9);
theta_OT = theta_OT.*(1 - abs(theta_OT)/pi*4/9);
theta_DOR = theta_DOR.*(1 - abs(theta_DOR)/pi*4/9);
%% saperation of data according to different v
v_list = [0.02, 0.11, 0.2];
data_processed.phi = {};
data_processed.lam = {};
data_processed.theta = {};
data_processed.u = {};
data_processed.v = {};
theta_predicted.CL = {};
theta_predicted.AP = {};
theta_predicted.PL = {};
theta_predicted.DO = {};
theta_predicted.OT = {};
theta_predicted.DOR = {};
for predator_speed = v_list
    data_processed.phi{end+1} = phi(v == predator_speed);
    data_processed.lam{end+1} = lam(v == predator_speed);
    data_processed.theta{end+1} = theta(v == predator_speed);
    data_processed.u{end+1} = u(v == predator_speed);
    data_processed.v{end+1} = v(v == predator_speed);
    theta_predicted.CL{end+1} = theta_CL(v == predator_speed);
    theta_predicted.AP{end+1} = theta_AP(v == predator_speed);
    theta_predicted.PL{end+1} = theta_PL(v == predator_speed);
    theta_predicted.DO{end+1} = theta_DO(v == predator_speed);
    theta_predicted.OT{end+1} = theta_OT(v == predator_speed);
    theta_predicted.DOR{end+1} = theta_DOR(v == predator_speed);
end
data_processed.phi{end+1} = phi;
data_processed.lam{end+1} = lam;
data_processed.theta{end+1} = theta;
data_processed.u{end+1} = u;
data_processed.v{end+1} = v;
theta_predicted.CL{end+1} = theta_CL;
theta_predicted.AP{end+1} = theta_AP;
theta_predicted.PL{end+1} = theta_PL;
theta_predicted.DO{end+1} = theta_DO;
theta_predicted.OT{end+1} = theta_OT;
theta_predicted.DOR{end+1} = theta_DOR;
% %%
% phi_02 = phi(v == 0.02);
% phi_11 = phi(v == 0.11);
% phi_20 = phi(v == 0.2);
% lam_02 = lam(v == 0.02);
% lam_11 = lam(v == 0.11);
% lam_20 = lam(v == 0.2);
% u_02 = u(v == 0.02);
% u_11 = u(v == 0.11);
% u_20 = u(v == 0.2);
% v_02 = v(v == 0.02);
% v_11 = v(v == 0.11);
% v_20 = v(v == 0.2);
% theta_02 = theta(v == 0.02);
% theta_11 = theta(v == 0.11);
% theta_20 = theta(v == 0.2);
%
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
% theta_DOR_02 = theta_DOR(v == 0.02);
% theta_DOR_11 = theta_DOR(v == 0.11);
% theta_DOR_20 = theta_DOR(v == 0.2);
%% histogram plots of each strategy VS experiment
fig = figure();
fig.Color = 'w';
fig.Position = [790 616 1075 211];
bin_edges = -pi-pi/36:pi/18:pi;
ax_CL = subplot(1,6,1,polaraxes);
ph_CL = polarhistogram(theta_CL,36,"BinEdges",bin_edges);

ax_AP = subplot(1,6,2,polaraxes);
ph_AP = polarhistogram(theta_AP,36,"BinEdges",bin_edges);

ax_PL = subplot(1,6,3,polaraxes);
ph_PL = polarhistogram(theta_PL,36,"BinEdges",bin_edges);

ax_DO = subplot(1,6,4,polaraxes);
ph_DO = polarhistogram(theta_DO,36,"BinEdges",bin_edges);

ax_OT = subplot(1,6,5,polaraxes);
ph_OT = polarhistogram(theta_OT,36,"BinEdges",bin_edges);

ax_exp = subplot(1,6,6,polaraxes);
ph_exp = polarhistogram(theta,36,"BinEdges",bin_edges);

ax_CL.FontSize = 8;
ax_AP.FontSize = 8;
ax_PL.FontSize = 8;
ax_DO.FontSize = 8;
ax_OT.FontSize = 8;
ax_exp.FontSize = 8;
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
theta_ticks = [0 45 90 135 180 225 270 315];
ax_CL.ThetaTick = theta_ticks;
ax_AP.ThetaTick = theta_ticks;
ax_PL.ThetaTick = theta_ticks;
ax_DO.ThetaTick = theta_ticks;
ax_OT.ThetaTick = theta_ticks;
ax_exp.ThetaTick = theta_ticks;
% show rad instead of degree
% theta_ticklabels = {'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
% ax_CL.ThetaTickLabel = theta_ticklabels;
% ax_AP.ThetaTickLabel = theta_ticklabels;
% ax_PL.ThetaTickLabel = theta_ticklabels;
% ax_DO.ThetaTickLabel = theta_ticklabels;
% ax_OT.ThetaTickLabel = theta_ticklabels;
% ax_exp.ThetaTickLabel = theta_ticklabels;
face_color = [0.2, 0.2, 0.2];
ph_CL.FaceColor = face_color;
ph_AP.FaceColor = face_color;
ph_PL.FaceColor = face_color;
ph_DO.FaceColor = face_color;
ph_OT.FaceColor = face_color;
ph_exp.FaceColor = face_color;

%% 3-D KL-Divergence per v
k = 5;
kl_reimagined.CL = zeros(1,4);
kl_reimagined.AP = zeros(1,4);
kl_reimagined.PL = zeros(1,4);
kl_reimagined.DO = zeros(1,4);
kl_reimagined.OT = zeros(1,4);
% kl_reimagined.DOR = zeros(1,4);
for j = 1:4
    X = [data_processed.theta{j}; data_processed.phi{j}; data_processed.lam{j}];
    Y = [theta_predicted.CL{j};  data_processed.phi{j}; data_processed.lam{j}];
    kl_reimagined.CL(j) = KLDivergence_angle(X,Y,k);
    Y(1,:) = theta_predicted.AP{j};
    kl_reimagined.AP(j) = KLDivergence_angle(X,Y,k);
    Y(1,:) = theta_predicted.PL{j};
    kl_reimagined.PL(j) = KLDivergence_angle(X,Y,k);
    Y(1,:) = theta_predicted.DO{j};
    kl_reimagined.DO(j) = KLDivergence_angle(X,Y,k);
    Y(1,:) = theta_predicted.OT{j};
    kl_reimagined.OT(j) = KLDivergence_angle(X,Y,k);
    Y(1,:) = theta_predicted.DOR{j};
    kl_reimagined.DOR(j) = KLDivergence_angle(X,Y,k);
end

% load('KLDivResults.mat')
K = [kl_reimagined.CL;kl_reimagined.AP;kl_reimagined.PL;kl_reimagined.DO;kl_reimagined.OT]';
fKL = figure;bp = bar(K); fKL.Position = [320 274 824 406];box off;
for j = 1:5
    bp(j).FaceColor = colorcode_strategy(j,:);bp(j).EdgeColor = 'none';
end
ax = gca; grid on; ax.GridLineStyle = '--';
ylim([-0.5,0.5]);
legend({'Contralateral','Antipodal','Parallel','Distance-optimal','Orthogonal'},'Location','bestoutside');
xticks('');
yticks(-0.5:0.5:2);
% title(['k = ' num2str(k)]);
%% inverse strategy --> finding best chi
k = 5;
chi = (0:36)*pi/72;
kl = zeros(4,length(chi));
for i = 1:4
    X = [data_processed.theta{i}; data_processed.phi{i}; data_processed.lam{i}];
    Y = X;
    branch = sign(data_processed.lam{i});
    for j = 1:length(chi)
        theta_DO = data_processed.phi{i} + data_processed.lam{i} + pi - branch*chi(j);
        theta_DO= wrapToPi(theta_DO);
        Y(1,:) = theta_DO;
        kl(i,j) = KLDivergence_angle(X,Y,k);
    end
end
[kl_min, index_min] = min(kl,[],2);
figure('Position',[921 400 783 472]);hold on;
lines = plot(chi,kl,'o-');
points = plot(chi(index_min),kl_min,'ks');
xlabel('$\chi$ [rad]');
xticks([0,pi/4,pi/2]); xticklabels({'0','$\pi/4$','$\pi/2$'});
ylabel('KL divergence');
title(['k = ' num2str(k)]);
legend('all','slow predator','medium speed predator','fast predator');

%% probabilistic evaluation

%% Setup Data
x_sample = cell(1,4);
Ns = zeros(1,4);
for i = 1:4
    Ns(i) = length(data_processed.phi{i});
    x_sample{i} = [data_processed.phi{i};...
        data_processed.lam{i};...
        data_processed.u{i};...
        data_processed.v{i};...
        data_processed.theta{i}];
end

%% Optimization option parameters
% Contralateral Strategy 
eta0_CL = [0.15*pi;...
    0.25*pi];

A_CL = [-1 0;...
    0 -1];

b_CL = [0;...
    0];
% Antipodal Strategy
eta0_AP = [0.22*pi;...  % sig_phi
    0.22*pi];    % sig_theta


A_AP = [-1  0 ;...
    0 -1 ];

b_AP = [0;...
    0];

% Parallel/Distance-optimal/Orthogonal Strategy
eta0_DO = [0.473078141407163;...  % sig_phi
    0.665034687972989;...  % sig_lam
    0.473088884833394];    % sig_theta


A_DO = [-1  0  0;
    0  -1  0;
    0   0 -1];

b_DO = [0;
    0;
    0];
%% Run Optimization

options = optimoptions('fmincon','Display','notify-detailed','MaxFunctionEvaluations',1e4);
%% Without bootstrapping
load('opt_initialize_reimagined.mat')
% eta_CL = [1;1];
% eta_AP = [1;1];
% eta_PL = [1;1;1];
% eta_DO = [1;1;1];
% eta_OT = [1;1;1];
eta_optimal = cell(1,4);
nll_optimal = cell(1,4);
for i =1:3
    i
    [eta_CL,nll_CL] = fmincon(@(eta)NLL_CL_reimagined(eta,x_sample{4}),eta0_CL,[],[],[],[],[0;0],[pi;pi],[],options);
    [eta_AP,nll_AP] = fmincon(@(eta)NLL_AP_reimagined(eta,x_sample{4}),eta0_AP,[],[],[],[],[0;0],[pi;pi],[],options);
    [eta_PL,nll_PL] = fmincon(@(eta)NLL_DO_reimagined(eta,x_sample{4},0),eta0_PL,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);
    [eta_DO,nll_DO] = fmincon(@(eta)NLL_DO_reimagined(eta,x_sample{4},acos(0.5)),eta0_DO,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);
    [eta_OT,nll_OT] = fmincon(@(eta)NLL_DO_reimagined(eta,x_sample{4},pi/2),eta0_OT,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);
    eta0_CL = eta_CL
    eta0_AP = eta_AP
    eta0_PL = eta_PL
    eta0_DO = eta_DO
    eta0_OT = eta_OT
    [eta_CL_02,nll_CL_02] = fmincon(@(eta)NLL_CL_reimagined(eta,x_sample{1}),eta0_CL_02,[],[],[],[],[0;0],[pi;pi],[],options);
    [eta_AP_02,nll_AP_02] = fmincon(@(eta)NLL_AP_reimagined(eta,x_sample{1}),eta0_AP_02,[],[],[],[],[0;0],[pi;pi],[],options);
    [eta_PL_02,nll_PL_02] = fmincon(@(eta)NLL_DO_reimagined(eta,x_sample{1},0),eta0_PL_02,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);
    [eta_DO_02,nll_DO_02] = fmincon(@(eta)NLL_DO_reimagined(eta,x_sample{1},acos(0.5)),eta0_DO_02,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);
    [eta_OT_02,nll_OT_02] = fmincon(@(eta)NLL_DO_reimagined(eta,x_sample{1},pi/2),eta0_OT_02,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);
    eta0_CL_02 = eta_CL_02
    eta0_AP_02 = eta_AP_02
    eta0_PL_02 = eta_PL_02
    eta0_DO_02 = eta_DO_02
    eta0_OT_02 = eta_OT_02
    [eta_CL_11,nll_CL_11] = fmincon(@(eta)NLL_CL_reimagined(eta,x_sample{2}),eta0_CL_11,[],[],[],[],[0;0],[pi;pi],[],options);
    [eta_AP_11,nll_AP_11] = fmincon(@(eta)NLL_AP_reimagined(eta,x_sample{2}),eta0_AP_11,[],[],[],[],[0;0],[pi;pi],[],options);
    [eta_PL_11,nll_PL_11] = fmincon(@(eta)NLL_DO_reimagined(eta,x_sample{2},0),eta0_PL_11,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);
    [eta_DO_11,nll_DO_11] = fmincon(@(eta)NLL_DO_reimagined(eta,x_sample{2},acos(0.5)),eta0_DO_11,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);
    [eta_OT_11,nll_OT_11] = fmincon(@(eta)NLL_DO_reimagined(eta,x_sample{2},pi/2),eta0_OT_11,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);
    eta0_CL_11 = eta_CL_11
    eta0_AP_11 = eta_AP_11
    eta0_PL_11 = eta_PL_11
    eta0_DO_11 = eta_DO_11
    eta0_OT_11 = eta_OT_11
    [eta_CL_20,nll_CL_20] = fmincon(@(eta)NLL_CL_reimagined(eta,x_sample{3}),eta0_CL_20,[],[],[],[],[0;0],[pi;pi],[],options);
    [eta_AP_20,nll_AP_20] = fmincon(@(eta)NLL_AP_reimagined(eta,x_sample{3}),eta0_AP_20,[],[],[],[],[0;0],[pi;pi],[],options);
    [eta_PL_20,nll_PL_20] = fmincon(@(eta)NLL_DO_reimagined(eta,x_sample{3},0),eta0_PL_20,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);
    [eta_DO_20,nll_DO_20] = fmincon(@(eta)NLL_DO_reimagined(eta,x_sample{3},acos(0.5)),eta0_DO_20,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);
    [eta_OT_20,nll_OT_20] = fmincon(@(eta)NLL_DO_reimagined(eta,x_sample{3},pi/2),eta0_OT_20,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);
    eta0_CL_20 = eta_CL_20
    eta0_AP_20 = eta_AP_20
    eta0_PL_20 = eta_PL_20
    eta0_DO_20 = eta_DO_20
    eta0_OT_20 = eta_OT_20
end
eta_optimal{j}.CL = eta_CL;
eta_optimal{j}.AP = eta_AP;
eta_optimal{j}.PL = eta_PL;
eta_optimal{j}.DO = eta_DO;
eta_optimal{j}.OT = eta_OT;
nll_optimal{j}.CL = nll_CL;
nll_optimal{j}.AP = nll_AP;
nll_optimal{j}.PL = nll_PL;
nll_optimal{j}.DO = nll_DO;
nll_optimal{j}.OT = nll_OT;
%% bootstrapping antipodal
load('opt_initialize_reimagined.mat');
load('booststrap_seed.mat')
%%
BSNum = 1000;
etaBS_AP_02  = zeros(2,BSNum);
etaBS_AP_11  = zeros(2,BSNum);
etaBS_AP_20  = zeros(2,BSNum);
etaBS_AP  = zeros(2,BSNum);
nllBS_AP_02 = zeros(1,BSNum);
nllBS_AP_11 = zeros(1,BSNum);
nllBS_AP_20 = zeros(1,BSNum);
nllBS_AP = zeros(1,BSNum);
for i = 1:BSNum
    %     r_02 = randi(Ns(1),[1,Ns(1)]);
    %     r_11 = randi(Ns(2),[1,Ns(2)]);
    %     r_20 = randi(Ns(3),[1,Ns(3)]);
    %     r = randi(Ns(4),[1,Ns(4)]);
    r_02 = seed_02(i,:);
    r_11 = seed_11(i,:);
    r_20 = seed_20(i,:);
    r = seed_combined(i,:);
    tic
    [etaBS_AP_02(:,i),nllBS_AP_02(i)] = fmincon(@(eta)NLL_AP_reimagined(eta,x_sample{1}(:,r_02)),eta0_AP_02,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);

    [etaBS_AP_11(:,i),nllBS_AP_11(i)] = fmincon(@(eta)NLL_AP_reimagined(eta,x_sample{2}(:,r_11)),eta0_AP_11,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);

    [etaBS_AP_20(:,i),nllBS_AP_20(i)] = fmincon(@(eta)NLL_AP_reimagined(eta,x_sample{3}(:,r_20)),eta0_AP_20,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);

    [etaBS_AP(:,i),nllBS_AP(i)] = fmincon(@(eta)NLL_AP_reimagined(eta,x_sample{4}(:,r)),eta0_AP,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);
    toc
    if exist('BSresults_reimagined_new.mat' ,'file')
        load('BSresults_reimagined_new.mat');
        eta_AP_02(:,end+1) = etaBS_AP_02(:,i);
        eta_AP_11(:,end+1) = etaBS_AP_11(:,i);
        eta_AP_20(:,end+1) = etaBS_AP_20(:,i);
        eta_AP(:,end+1) = etaBS_AP(:,i);
        nll_AP_02(end+1) = nllBS_AP_02(i);
        nll_AP_11(end+1) = nllBS_AP_11(i);
        nll_AP_20(end+1) = nllBS_AP_20(i);
        nll_AP(end+1) = nllBS_AP(i);
        save BSresults_reimagined_new.mat eta_AP_02 eta_AP_11 eta_AP_20 eta_AP nll_AP_02 nll_AP_11 nll_AP_20 nll_AP -append
    else
        eta_AP_02 = etaBS_AP_02(:,i);
        eta_AP_11 = etaBS_AP_11(:,i);
        eta_AP_20 = etaBS_AP_20(:,i);
        eta_AP = etaBS_AP(:,1);
        nll_AP_02 = nllBS_AP_02(i);
        nll_AP_11 = nllBS_AP_11(i);
        nll_AP_20 = nllBS_AP_20(i);
        nll_AP = nllBS_AP(i);
        save BSresults_reimagined_new.mat eta_AP_02 eta_AP_11 eta_AP_20 eta_AP nll_AP_02 nll_AP_11 nll_AP_20 nll_AP
    end
    
end
%% bootstrapping Orthogonal
load('opt_initialize_reimagined.mat');

chi = pi/2;
BSNum = 73;
etaBS_OT_02  = zeros(3,BSNum);
etaBS_OT_11  = zeros(3,BSNum);
etaBS_OT_20  = zeros(3,BSNum);
etaBS_OT  = zeros(3,BSNum);
nllBS_OT_02 = zeros(1,BSNum);
nllBS_OT_11 = zeros(1,BSNum);
nllBS_OT_20 = zeros(1,BSNum);
nllBS_OT = zeros(1,BSNum);
% eta_OT_02=[];eta_OT_11=[];eta_OT_20=[];eta_OT=[];
% nll_OT_02=[];nll_OT_11=[];nll_OT_20=[];nll_OT=[];
% save BSresults_reimagined.mat eta_OT_02 eta_OT_11 eta_OT_20 eta_OT nll_OT_02 nll_OT_11 nll_OT_20 nll_OT -append
for i = 1:BSNum
    r_02 = seed_02(i,:);
    r_11 = seed_11(i,:);
    r_20 = seed_20(i,:);
    r = seed_combined(i,:);
    tic
    [etaBS_OT_02(:,i),nllBS_OT_02(i)] = fmincon(@(eta)NLL_DO_reimagined(eta,x_sample{1}(:,r_02),chi),eta0_OT_02,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);

    [etaBS_OT_11(:,i),nllBS_OT_11(i)] = fmincon(@(eta)NLL_DO_reimagined(eta,x_sample{2}(:,r_11),chi),eta0_OT_11,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);

    [etaBS_OT_20(:,i),nllBS_OT_20(i)] = fmincon(@(eta)NLL_DO_reimagined(eta,x_sample{3}(:,r_20),chi),eta0_OT_20,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);

    [etaBS_OT(:,i),nllBS_OT(i)] = fmincon(@(eta)NLL_DO_reimagined(eta,x_sample{4}(:,r),chi),eta0_OT,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);
    toc
    if exist('BSresults_reimagined_new_OT.mat','file')
        load('BSresults_reimagined_new_OT.mat');
        eta_OT_02(:,end+1) = etaBS_OT_02(:,i);
        eta_OT_11(:,end+1) = etaBS_OT_11(:,i);
        eta_OT_20(:,end+1) = etaBS_OT_20(:,i);
        eta_OT(:,end+1) = etaBS_OT(:,i);
        nll_OT_02(end+1) = nllBS_OT_02(i);
        nll_OT_11(end+1) = nllBS_OT_11(i);
        nll_OT_20(end+1) = nllBS_OT_20(i);
        nll_OT(end+1) = nllBS_OT(i);
        save BSresults_reimagined_new_OT.mat eta_OT_02 eta_OT_11 eta_OT_20 eta_OT nll_OT_02 nll_OT_11 nll_OT_20 nll_OT -append
    else
        eta_OT_02 = etaBS_OT_02(:,i);
        eta_OT_11 = etaBS_OT_11(:,i);
        eta_OT_20 = etaBS_OT_20(:,i);
        eta_OT = etaBS_OT(:,i);
        nll_OT_02 = nllBS_OT_02(i);
        nll_OT_11 = nllBS_OT_11(i);
        nll_OT_20 = nllBS_OT_20(i);
        nll_OT = nllBS_OT(i);
        save BSresults_reimagined_new_OT.mat eta_OT_02 eta_OT_11 eta_OT_20 eta_OT nll_OT_02 nll_OT_11 nll_OT_20 nll_OT
    end
    
end

%% bootstrapping Distance-optimal with u/v=0.5
load('opt_initialize_reimagined.mat');
chi = acos(0.5);
BSNum = 19;
% etaBS_DO_02  = zeros(3,BSNum);
% etaBS_DO_11  = zeros(3,BSNum);
% etaBS_DO_20  = zeros(3,BSNum);
% etaBS_DO  = zeros(3,BSNum);
% nllBS_DO_02 = zeros(1,BSNum);
% nllBS_DO_11 = zeros(1,BSNum);
% nllBS_DO_20 = zeros(1,BSNum);
% nllBS_DO = zeros(1,BSNum);
eta_DO_02=[];eta_DO_11=[];eta_DO_20=[];eta_DO=[];
nll_DO_02=[];nll_DO_11=[];nll_DO_20=[];nll_DO=[];
% save BSresults_reimagined.mat eta_DO_02 eta_DO_11 eta_DO_20 eta_DO nll_DO_02 nll_DO_11 nll_DO_20 nll_DO -append
% parpool('threads');
for i = 1:BSNum
    r_02 = seed_02(i,:);
    r_11 = seed_11(i,:);
    r_20 = seed_20(i,:);
    r = seed_combined(i,:);
    tic
    [etaBS_DO_02,nllBS_DO_02] = fmincon(@(eta)NLL_DO_reimagined(eta,x_02(:,r_02),chi),eta0_DO_02,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);

    [etaBS_DO_11,nllBS_DO_11] = fmincon(@(eta)NLL_DO_reimagined(eta,x_11(:,r_11),chi),eta0_DO_11,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);

    [etaBS_DO_20,nllBS_DO_20] = fmincon(@(eta)NLL_DO_reimagined(eta,x_20(:,r_20),chi),eta0_DO_20,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);

    [etaBS_DO,nllBS_DO] = fmincon(@(eta)NLL_DO_reimagined(eta,x(:,r),chi),eta0_DO,[],[],[],[],[1e-2;1e-2;1e-2],[pi;pi;pi],[],options);
    toc
    load('BSresults_reimagined_new.mat');
    eta_DO_02(:,end+1) = etaBS_DO_02;
    eta_DO_11(:,end+1) = etaBS_DO_11;
    eta_DO_20(:,end+1) = etaBS_DO_20;
    eta_DO(:,end+1) = etaBS_DO;
    nll_DO_02(end+1) = nllBS_DO_02;
    nll_DO_11(end+1) = nllBS_DO_11;
    nll_DO_20(end+1) = nllBS_DO_20;
    nll_DO(end+1) = nllBS_DO;
    save BSresults_reimagined_new.mat eta_DO_02 eta_DO_11 eta_DO_20 eta_DO nll_DO_02 nll_DO_11 nll_DO_20 nll_DO -append
end
%% Compute AIC and Make the figure

load('BSresultsAll.mat');
load('BSresults_reimagined.mat');

N = length(nllopt_DO);
mAIC_02 = zeros(6,N);
mAIC_11 = zeros(6,N);
mAIC_20 = zeros(6,N);
mAIC_All = zeros(6,N);
sAIC_02 = zeros(6,N);
sAIC_11 = zeros(6,N);
sAIC_20 = zeros(6,N);
sAIC_All = zeros(6,N);

Datasize = [251;233;215;699];
AICopt_DO = nllopt_DO.*Datasize*2 + 4;
AICopt_OT = nllopt_OT.*Datasize*2 + 4;
AICopt_AP = nllopt_AP.*Datasize*2 + 2;

AICopt_DO_reimagined(1,:) = nll_DO_02*Datasize(1)*2 + 4;
AICopt_DO_reimagined(2,:) = nll_DO_11*Datasize(2)*2 + 4;
AICopt_DO_reimagined(3,:) = nll_DO_20*Datasize(3)*2 + 4;
AICopt_DO_reimagined(4,:) = nll_DO*Datasize(4)*2 + 4;

AICopt_OT_reimagined(1,:) = nll_OT_02*Datasize(1)*2 + 4;
AICopt_OT_reimagined(2,:) = nll_OT_11*Datasize(2)*2 + 4;
AICopt_OT_reimagined(3,:) = nll_OT_20*Datasize(3)*2 + 4;
AICopt_OT_reimagined(4,:) = nll_OT*Datasize(4)*2 + 4;

AICopt_AP_reimagined(1,:) = nll_AP_02*Datasize(1)*2 + 2;
AICopt_AP_reimagined(2,:) = nll_AP_11*Datasize(2)*2 + 2;
AICopt_AP_reimagined(3,:) = nll_AP_20*Datasize(3)*2 + 2;
AICopt_AP_reimagined(4,:) = nll_AP*Datasize(4)*2 + 2;

for i = 1:length(AICopt_DO)
    mAIC_02(1,i) = mean(AICopt_DO(1,1:i)); % Distance-optimal
    mAIC_11(1,i) = mean(AICopt_DO(2,1:i)); % Distance-optimal
    mAIC_20(1,i) = mean(AICopt_DO(3,1:i)); % Distance-optimal
    mAIC_All(1,i) = mean(AICopt_DO(4,1:i)); % Distance-optimal
    sAIC_02(1,i) = std(AICopt_DO(1,1:i)); % Distance-optimal
    sAIC_11(1,i) = std(AICopt_DO(2,1:i)); % Distance-optimal
    sAIC_20(1,i) = std(AICopt_DO(3,1:i)); % Distance-optimal
    sAIC_All(1,i) = std(AICopt_DO(4,1:i)); % Distance-optimal
end
for i = 1:length(AICopt_DO_reimagined)
    mAIC_02(2,i) = mean(AICopt_DO_reimagined(1,1:i));
    mAIC_11(2,i) = mean(AICopt_DO_reimagined(2,1:i));
    mAIC_20(2,i) = mean(AICopt_DO_reimagined(3,1:i));
    mAIC_All(2,i) = mean(AICopt_DO_reimagined(4,1:i));
    sAIC_02(2,i) = std(AICopt_DO_reimagined(1,1:i));
    sAIC_11(2,i) = std(AICopt_DO_reimagined(2,1:i));
    sAIC_20(2,i) = std(AICopt_DO_reimagined(3,1:i));
    sAIC_All(2,i) = std(AICopt_DO_reimagined(4,1:i));
end
for i = 1:length(AICopt_OT)
    mAIC_02(3,i) = mean(AICopt_OT(1,1:i)); % Orthogonal
    mAIC_11(3,i) = mean(AICopt_OT(2,1:i)); % Orthogonal
    mAIC_20(3,i) = mean(AICopt_OT(3,1:i)); % Orthogonal
    mAIC_All(3,i) = mean(AICopt_OT(4,1:i)); % Orthogonal
    sAIC_02(3,i) = std(AICopt_OT(1,1:i)); % Orthogonal
    sAIC_11(3,i) = std(AICopt_OT(2,1:i)); % Orthogonal
    sAIC_20(3,i) = std(AICopt_OT(3,1:i)); % Orthogonal
    sAIC_All(3,i) = std(AICopt_OT(4,1:i)); % Orthogonal
end
for i = 1:length(AICopt_OT_reimagined)
    mAIC_02(4,i) = mean(AICopt_OT_reimagined(1,1:i));
    mAIC_11(4,i) = mean(AICopt_OT_reimagined(2,1:i));
    mAIC_20(4,i) = mean(AICopt_OT_reimagined(3,1:i));
    mAIC_All(4,i) = mean(AICopt_OT_reimagined(4,1:i));
    sAIC_02(4,i) = std(AICopt_OT_reimagined(1,1:i));
    sAIC_11(4,i) = std(AICopt_OT_reimagined(2,1:i));
    sAIC_20(4,i) = std(AICopt_OT_reimagined(3,1:i));
    sAIC_All(4,i) = std(AICopt_OT_reimagined(4,1:i));
end
for i = 1:length(AICopt_AP)
    mAIC_02(5,i) = mean(AICopt_AP(1,1:i)); % Antipodal
    mAIC_11(5,i) = mean(AICopt_AP(2,1:i)); % 
    mAIC_20(5,i) = mean(AICopt_AP(3,1:i)); % 
    mAIC_All(5,i) = mean(AICopt_AP(4,1:i)); %
    sAIC_02(5,i) = std(AICopt_AP(1,1:i)); % 
    sAIC_11(5,i) = std(AICopt_AP(2,1:i)); % 
    sAIC_20(5,i) = std(AICopt_AP(3,1:i)); % 
    sAIC_All(5,i) = std(AICopt_AP(4,1:i)); % 
end
for i = 1:length(AICopt_AP_reimagined)
    mAIC_02(6,i) = mean(AICopt_AP_reimagined(1,1:i));
    mAIC_11(6,i) = mean(AICopt_AP_reimagined(2,1:i));
    mAIC_20(6,i) = mean(AICopt_AP_reimagined(3,1:i));
    mAIC_All(6,i) = mean(AICopt_AP_reimagined(4,1:i));
    sAIC_02(6,i) = std(AICopt_AP_reimagined(1,1:i));
    sAIC_11(6,i) = std(AICopt_AP_reimagined(2,1:i));
    sAIC_20(6,i) = std(AICopt_AP_reimagined(3,1:i));
    sAIC_All(6,i) = std(AICopt_AP_reimagined(4,1:i));
end

AICm_norm_02_DO = nonzeros(mAIC_02(1,:))/Datasize(1);
AICm_norm_02_DO_reimagined = nonzeros(mAIC_02(2,:))/Datasize(1);
AICm_norm_02_OT = nonzeros(mAIC_02(3,:))/Datasize(1);
AICm_norm_02_OT_reimagined = nonzeros(mAIC_02(4,:))/Datasize(1);
AICm_norm_02_AP = nonzeros(mAIC_02(5,:))/Datasize(1);
AICm_norm_02_AP_reimagined = nonzeros(mAIC_02(6,:))/Datasize(1);

AICm_norm_11_DO = nonzeros(mAIC_11(1,:))/Datasize(2);
AICm_norm_11_DO_reimagined = nonzeros(mAIC_11(2,:))/Datasize(2);
AICm_norm_11_OT = nonzeros(mAIC_11(3,:))/Datasize(2);
AICm_norm_11_OT_reimagined = nonzeros(mAIC_11(4,:))/Datasize(2);
AICm_norm_11_AP = nonzeros(mAIC_11(5,:))/Datasize(2);
AICm_norm_11_AP_reimagined = nonzeros(mAIC_11(6,:))/Datasize(2);

AICm_norm_20_DO = nonzeros(mAIC_20(1,:))/Datasize(3);
AICm_norm_20_DO_reimagined = nonzeros(mAIC_20(2,:))/Datasize(3);
AICm_norm_20_OT = nonzeros(mAIC_20(3,:))/Datasize(3);
AICm_norm_20_OT_reimagined = nonzeros(mAIC_20(4,:))/Datasize(3);
AICm_norm_20_AP = nonzeros(mAIC_20(5,:))/Datasize(3);
AICm_norm_20_AP_reimagined = nonzeros(mAIC_20(6,:))/Datasize(3);

AICm_norm_All_DO = nonzeros(mAIC_All(1,:))/Datasize(4);
AICm_norm_All_DO_reimagined = nonzeros(mAIC_All(2,:))/Datasize(4);
AICm_norm_All_OT = nonzeros(mAIC_All(3,:))/Datasize(4);
AICm_norm_All_OT_reimagined = nonzeros(mAIC_All(4,:))/Datasize(4);
AICm_norm_All_AP = nonzeros(mAIC_All(5,:))/Datasize(4);
AICm_norm_All_AP_reimagined = nonzeros(mAIC_All(6,:))/Datasize(4);

AICs_norm_02_DO = nonzeros(sAIC_02(1,:))/Datasize(1);
AICs_norm_02_DO_reimagined = nonzeros(sAIC_02(2,:))/Datasize(1);
AICs_norm_02_OT = nonzeros(sAIC_02(3,:))/Datasize(1);
AICs_norm_02_OT_reimagined = nonzeros(sAIC_02(4,:))/Datasize(1);
AICs_norm_02_AP = nonzeros(sAIC_02(5,:))/Datasize(1);
AICs_norm_02_AP_reimagined = nonzeros(sAIC_02(6,:))/Datasize(1);

AICs_norm_11_DO = nonzeros(sAIC_11(1,:))/Datasize(2);
AICs_norm_11_DO_reimagined = nonzeros(sAIC_11(2,:))/Datasize(2);
AICs_norm_11_OT = nonzeros(sAIC_11(3,:))/Datasize(2);
AICs_norm_11_OT_reimagined = nonzeros(sAIC_11(4,:))/Datasize(2);
AICs_norm_11_AP = nonzeros(sAIC_11(5,:))/Datasize(2);
AICs_norm_11_AP_reimagined = nonzeros(sAIC_11(6,:))/Datasize(2);

AICs_norm_20_DO = nonzeros(sAIC_20(1,:))/Datasize(3);
AICs_norm_20_DO_reimagined = nonzeros(sAIC_20(2,:))/Datasize(3);
AICs_norm_20_OT = nonzeros(sAIC_20(3,:))/Datasize(3);
AICs_norm_20_OT_reimagined = nonzeros(sAIC_20(4,:))/Datasize(3);
AICs_norm_20_AP = nonzeros(sAIC_20(5,:))/Datasize(3);
AICs_norm_20_AP_reimagined = nonzeros(sAIC_20(6,:))/Datasize(3);

AICs_norm_All_DO = nonzeros(sAIC_All(1,:))/Datasize(4);
AICs_norm_All_DO_reimagined = nonzeros(sAIC_All(2,:))/Datasize(4);
AICs_norm_All_OT = nonzeros(sAIC_All(3,:))/Datasize(4);
AICs_norm_All_OT_reimagined = nonzeros(sAIC_All(4,:))/Datasize(4);
AICs_norm_All_AP = nonzeros(sAIC_All(5,:))/Datasize(4);
AICs_norm_All_AP_reimagined = nonzeros(sAIC_All(6,:))/Datasize(4);

AICm_norm_02 = [AICm_norm_02_DO(end),AICm_norm_02_DO_reimagined(end),...
    AICm_norm_02_OT(end),AICm_norm_02_OT_reimagined(end),...
    AICm_norm_02_AP(end),AICm_norm_02_AP_reimagined(end)];
AICm_norm_11 = [AICm_norm_11_DO(end),AICm_norm_11_DO_reimagined(end),...
    AICm_norm_11_OT(end),AICm_norm_11_OT_reimagined(end),...
    AICm_norm_11_AP(end),AICm_norm_11_AP_reimagined(end)];
AICm_norm_20 = [AICm_norm_20_DO(end),AICm_norm_20_DO_reimagined(end),...
    AICm_norm_20_OT(end),AICm_norm_20_OT_reimagined(end),...
    AICm_norm_11_AP(end),AICm_norm_11_AP_reimagined(end)];
AICm_norm_All = [AICm_norm_All_DO(end),AICm_norm_All_DO_reimagined(end),...
    AICm_norm_All_OT(end),AICm_norm_All_OT_reimagined(end),...
    AICm_norm_All_AP(end),AICm_norm_All_AP_reimagined(end)];
AICs_norm_02 = [AICs_norm_02_DO(end),AICs_norm_02_DO_reimagined(end),...
    AICs_norm_02_OT(end),AICs_norm_02_OT_reimagined(end),...
    AICs_norm_02_AP(end),AICs_norm_02_AP_reimagined(end)];
AICs_norm_11 = [AICs_norm_11_DO(end),AICs_norm_11_DO_reimagined(end),...
    AICs_norm_11_OT(end),AICs_norm_11_OT_reimagined(end),...
    AICs_norm_11_AP(end),AICs_norm_11_AP_reimagined(end)];
AICs_norm_20 = [AICs_norm_20_DO(end),AICs_norm_20_DO_reimagined(end),...
    AICs_norm_20_OT(end),AICs_norm_20_OT_reimagined(end),...
    AICs_norm_20_AP(end),AICs_norm_20_AP_reimagined(end)];
AICs_norm_All = [AICs_norm_All_DO(end),AICs_norm_All_DO_reimagined(end),...
    AICs_norm_All_OT(end),AICs_norm_All_OT_reimagined(end),...
    AICs_norm_All_AP(end),AICs_norm_All_AP_reimagined(end)];

x = [1:6 , 8:13, 15:20, 22:27];
meanArray = [AICm_norm_02-min(AICm_norm_02), AICm_norm_11-min(AICm_norm_11),...
    AICm_norm_20-min(AICm_norm_20), AICm_norm_All-min(AICm_norm_All)];
stdArray = [AICs_norm_02, AICs_norm_11, AICs_norm_20, AICs_norm_All];

figure('position',[253 1177 321 182]);
eb = errorbar(x, meanArray,stdArray,'marker','s','linestyle','none','markersize',4,'capsize',0);
box off;
xticks(3.5:7:27);
xlabel('dataset');
ylim([-0.2,0.5]);
ylabel('$\Delta$AIC$/N$');
yticks(-0.2:0.1:0.5);