clear
close all
clc

%% Load Data

load('evasionData.mat')

%% Sort Data

phi = data.Phi(isfinite(data.U));
u   = data.U(isfinite(data.U));
v   = data.V(isfinite(data.U))*0.01;
lam = data.Lambda(isfinite(data.U));
theta = data.Delta(isfinite(data.U));

%% favorable sample for S2

% phi_2 = phi(abs(cos(phi - theta + pi)-1)<0.00001);
% lam_2 = lam(abs(cos(phi - theta + pi)-1)<0.00001);
% theta_2 = theta(abs(cos(phi - theta + pi)-1)<0.00001);
% 
% phi = phi_2;
% lam = lam_2;
% theta = theta_2;
%% favorable sample for S3

phi_select = phi(abs(cos(phi + lam - theta + pi - sign(lam)*pi/2)-1)<0.000005);
lam_select = lam(abs(cos(phi + lam - theta + pi - sign(lam)*pi/2)-1)<0.000005);
theta_select = theta(abs(cos(phi + lam - theta + pi - sign(lam)*pi/2)-1)<0.000005);
%%

for k = 1:length(lam_select)
    figure,quiver(0,0,1,0,'LineWidth',2);
    hold on;
    quiver(2*cos(phi_select(k)), 2*sin(phi_select(k)), cos(phi_select(k)+lam_select(k)),sin(phi_select(k)+lam_select(k)),'LineWidth',2);
    quiver(0,0,cos(theta_select(k)),sin(theta_select(k)),'LineWidth',2,'LineStyle','--');
    xlim([-3,3]);
    ylim([-3,3]);
    axis equal;
end