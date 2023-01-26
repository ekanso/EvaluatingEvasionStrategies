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

%% swim circle
figure(), hold on;
load('swim1.mat');
% plot(var(:,1),var(:,2));
plot(rx_mean_history(1:180)-rx_mean_history(1),ry_mean_history(1:180)-ry_mean_history(1));
load('swim2.mat');
% plot(var(:,1),var(:,2));
plot(rx_mean_history(1:180)-rx_mean_history(1),ry_mean_history(1:180)-ry_mean_history(1));
load('swim3.mat');
% plot(var(:,1),var(:,2));
plot(rx_mean_history(1:180)-rx_mean_history(1),ry_mean_history(1:180)-ry_mean_history(1));
load('swim4.mat');
% plot(var(:,1),var(:,2));
plot(rx_mean_history(1:180)-rx_mean_history(1),ry_mean_history(1:180)-ry_mean_history(1));
axis equal; grid on;
xlim([-0.1,0.9]);
ylim([-0.5,0.5]);

%% Turn circle
figure(), hold on;
load('turn1.mat');
% plot(var(:,1),var(:,2));
plot(rx_mean_history(1:180)-rx_mean_history(1),ry_mean_history(1:180)-ry_mean_history(1));
load('turn2.mat');
% plot(var(:,1),var(:,2));
plot(rx_mean_history(1:180)-rx_mean_history(1),ry_mean_history(1:180)-ry_mean_history(1));
load('turn3.mat');
% plot(var(:,1),var(:,2));
plot(rx_mean_history(1:180)-rx_mean_history(1),ry_mean_history(1:180)-ry_mean_history(1));
load('turn4.mat');
% plot(var(:,1),var(:,2));
plot(rx_mean_history(1:180)-rx_mean_history(1),ry_mean_history(1:180)-ry_mean_history(1));
axis equal; grid on;
xlim([-0.2,0.6]);
ylim([-0.6,0.2]);

%% patched circles to turn 60 deg
figure(), hold on;
load('turnPlusSwim.mat');

plot(rx_mean_history(1:180)-rx_mean_history(1),ry_mean_history(1:180)-ry_mean_history(1));load('turn2.mat');
axis equal; grid off;box on;
xlim([-0.1,0.8]);
ylim([-0.1,0.2]);