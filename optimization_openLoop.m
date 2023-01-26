clear
close all
clc
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

Ns = length(phi);

x = [phi;
    lam;
    u;
    v;
    theta];
%%
% map = zeros(101);
% for i = 1:101
%     tic
%     for j = 1:101
%         map(i,j) = NLL_S3([pi/100*(i-1);etaopt(2);pi/100*(j-1)],x);
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
diary mapS3
% map_phi = zeros(11,51,51);
% map_lambda = zeros(11,51,51);
map = zeros(71,71,71);
for k = 1:71
    for i = 1:71
        tic
        %         for j = 1:51
        %             map_phi(k,i,j) = NLL_S3([pi*(i-1)/100; etaopt(2)*(4/5 + (k-1)/25); pi/100*(j-1)],x);
        %         end
        %
        %         for j = 1:51
        %             map_lambda(k,i,j) = NLL_S3([etaopt(1)*(4/5 + (k-1)/25); pi/100*(i-1); 0.1 + pi/100*(j-1)],x);
        %         end
        for j = 1:71
            map(i,j,k) = NLL_DO([pi*(i-1)/100; pi/100*(j-1); pi/100*(k-1)],x,pi/2);
        end
        toc
    end
end
% map_v = zeros(25);
% for i = 1:25
%     tic
%     for j = 1:25
%         map_v(i,j) = NLL_S4([0.5174; 1.2212; 0.15 + 0.6*(i-1)/24; 0.1 + 0.6/24*(j-1)],x);
%     end
%     toc
%     i
% end
% map_v
diary off
save s3_map.mat map
%%
cmap=hot(30);
[M,cf] = contourf(0:pi/120:pi/2,0:pi/120:pi,map(:,:,5).','ShowText','off','EdgeColor','none');
view(2);
cf.LevelList = [1.51:0.01:1.70 1.72:0.02:(-log(1/2/pi))];
cf.LineStyle = '--';
colormap(cmap);
colorbar('Ticks',[min(min(map(:,:,5))),1.9]);
caxis([min(min(map(:,:,5))) 1.9]);
title('Strategy 3', 'Interpreter','latex');
ax = gca;
ax.XTick = [0 pi/2 pi];
ax.YTick = [0 pi/2 pi];
ax.XTickLabel = {'0','$\pi/2$','$\pi$'};
ax.YTickLabel = {'0','$\pi/2$','$\pi$'};
ax.XLabel.String = '$\sigma_\Lambda$';
ax.XLabel.Interpreter = 'latex';
ax.FontSize = 18;
ax.TickLabelInterpreter = 'latex';
ax.YLabel.String = '$\sigma_\Theta$';
ax.YLabel.Interpreter = 'latex';
%%
N_gridLam = 121;
N_gridPhi = 61;
chi = (1:10)*pi/2/11;
N_slice = length(chi);
% map = zeros(N_gridPhi,N_gridLam,);
for j = 5:N_slice
for k = 1:N_gridLam
    tic
    for i = 1:N_gridPhi
        map(i,k,j) = NLL_DO([pi*(i-1)/N_gridPhi/2; pi*(k-1)/N_gridLam],x_20,chi(j));
    end
    toc
    save openloop_maps_Chi.mat map
end
end
%%
cmap=hot(30);
chi = (1:10)*pi/2/11;
for j = 1:size(map,3)
figure;
[M,cf] = contourf(0:pi/120:pi/2,0:pi/120:pi,map(:,:,j)','ShowText','off','EdgeColor','none');
view(2);
cf.LevelList = [cf.LevelList(1):0.01:(-log(1/2/pi))];
cf.LineStyle = '--';
colormap(cmap);
% colorbar('Ticks',[min(min(map(:,:,j))),-log(1/2/pi)]);
caxis([min(min(map(:,:,j))) -log(1/2/pi)]);
title(['Distance-optimal strategy with $\chi=$' num2str(chi(j))], 'Interpreter','latex');
ax = gca;
ax.XTick = [0 pi/2];
ax.YTick = [0 pi/2 pi];
ax.XTickLabel = {'0','$\pi/2$'};
ax.YTickLabel = {'0','$\pi/2$','$\pi$'};
ax.YLabel.String = '$\sigma_\Lambda$';
ax.XLabel.Interpreter = 'latex';
ax.FontSize = 18;
ax.TickLabelInterpreter = 'latex';
ax.XLabel.String = '$\sigma_\Phi, \sigma_\Theta$';
ax.YLabel.Interpreter = 'latex';
axis equal
end
%% 
L1 = log(pX_DO_mix_test(x_20,[0.48;0.54;0.48],pi/2));
sum(L1)
L2 = log(pX_DO_mix_test(x_20,[0.48;pi;0.48],pi/2));
L3 = log(pX_DO_mix_test(x_20,[pi;pi;pi],pi/2));
sum(L2)