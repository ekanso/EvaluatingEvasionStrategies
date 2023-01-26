% use beta_m to denote the orientation of the middle body, in order to
% avoid confliction with the built-in function beta()

clear;
close all;
delta_his = [];
delta_mean_his = [];
beta_his = [];
beta_mean_his = [];
ii = 1;
%%
frame_num = 121;
time_expect = 2*pi; % total simulation time
act_mode = 'CStartFit';
%% Visualization
fig = figure('Position',[225 377 1150 420]); hold all; grid off;
% fig.Color = 'none';
axis equal
xlim([-0.5,1.5]);
ylim([-0.5,2]);
%% description of the 3 links
a = 1;
b = 0.2;
c = 1;

L_c = 6*a;% to make the size dimensionless
a = a/L_c;
b = b/L_c;
c = c/L_c;

rho = 1000;
M_c = 4/3*rho*(a*b*c+a*b*c+a*b*c)*pi;

%% body mass and added mass
m = zeros(3,1);
J = zeros(2,1);

m(1) = a*b*c/(a*b*c + a*b*c + a*b*c);
J(1) = 1/5*(a^2+b^2)*m(1);

% calculate the added mass and moment of inertia
[m(2),m(3),J(2)] = getAddedMass(a,b,c);

%% initialization of the variations
[xi,Ixi] = initialConfig(act_mode,ii);

rx = Ixi(1);
ry = Ixi(2);
beta_m = Ixi(3);
alpha1 = Ixi(4);
alpha2 = Ixi(5);

r2x = rx - cos(beta_m)*a - cos(beta_m-alpha2)*a;
r2y = ry - sin(beta_m)*a - sin(beta_m-alpha2)*a;
r1x = rx + cos(beta_m)*a + cos(beta_m+alpha1)*a;
r1y = ry + sin(beta_m)*a + sin(beta_m+alpha1)*a;


xi_history = xi;
Ixi_history = Ixi;

beta_mean_ini = beta_m - alpha2/3 + alpha1/3;


%% Ignore body mass
m(1) = 0;
J(1) = 0;
%% time integration
var = Ixi(1:3);
Momt = [0;0;0];
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,var] = ode45(@(t,var) threelink_dynamics(t,act_mode,var,a,m,J,ii),linspace(0,time_expect,frame_num),var,opts);
alpha1_history = zeros(length(t),1);
alpha2_history = zeros(length(t),1);
for t_index = 1:length(t)

    [alpha1_history(t_index),alpha2_history(t_index),~,~] = prescribedAngle(t(t_index),act_mode,ii);
end


rx_history = var(:,1);
ry_history = var(:,2);
beta_m_history = var(:,3);
beta2_history = beta_m_history-alpha2_history;
beta1_history = beta_m_history+alpha1_history;
r2x_history = rx_history - cos(beta_m_history)*a - cos(beta2_history)*a;
r2y_history = ry_history - sin(beta_m_history)*a - sin(beta2_history)*a;
r1x_history = rx_history + cos(beta_m_history)*a + cos(beta1_history)*a;
r1y_history = ry_history + sin(beta_m_history)*a + sin(beta1_history)*a;
rx_mean_history = (rx_history+r2x_history+r1x_history)/3;
ry_mean_history = (ry_history+r2y_history+r1y_history)/3;

%% snapshots
for index = 1:100:frame_num
    Ixi = var(index,:)';
    [alpha1,alpha2,Dalpha1,Dalph2] = prescribedAngle(t(index),act_mode,ii);

    rx = Ixi(1);
    ry = Ixi(2);
    beta_m = Ixi(3);
    beta1 = beta_m+alpha1;
    beta2 = beta_m-alpha2;
    %     r2x = rx - cos(beta_m)*a - cos(beta2)*a;
    %     r2y = ry - sin(beta_m)*a - sin(beta2)*a;
    %     r1x = rx + cos(beta_m)*a + cos(beta1)*a;
    %     r1y = ry + sin(beta_m)*a + sin(beta1)*a;
    %     x_mean = (rx+r2x+r1x)/3;
    %     y_mean = (ry+r2y+r1y)/3;
    disp([rx*10/a/2, ry*10/a/2, beta_m*180/pi])
    disp(rad2deg(atan2(alpha2,alpha1)))
    make3linkfish(a,b,c,rx,ry,beta_m,alpha1,alpha2,'2d');
end
plot(rx_mean_history,ry_mean_history,'-')
%% video
v = VideoWriter(['fishmovie' act_mode '.avi']);
v.Quality = 100;
v.FrameRate = 30;
axis off
open(v)
Ixi = var(1,:)';
[alpha1,alpha2,Dalpha1,Dalph2] = prescribedAngle(t(1),act_mode,ii);

rx = Ixi(1);
ry = Ixi(2);
beta_m = Ixi(3);
beta1 = beta_m+alpha1;
beta2 = beta_m-alpha2;
h_fish = make3linkfish(a,b,c,rx,ry,beta_m,alpha1,alpha2,'2d');
traj = plot(rx_mean_history(1),ry_mean_history(1),'.r');
frame = getframe(gcf);
writeVideo(v,frame)
for index = 2:frame_num
    Ixi = var(index,:)';
    [alpha1,alpha2,Dalpha1,Dalph2] = prescribedAngle(t(index),act_mode,ii);
    rx = Ixi(1);
    ry = Ixi(2);
    beta_m = Ixi(3);
    beta1 = beta_m+alpha1;
    beta2 = beta_m-alpha2;


    h_fish = make3linkfish(a,b,c,rx,ry,beta_m,alpha1,alpha2,'2d',h_fish);
    traj.XData = rx_mean_history(1:index);
    traj.YData = ry_mean_history(1:index);
%     plot(rx_mean_history(index),ry_mean_history(index),'.r')
    frame = getframe(gcf);
    writeVideo(v,frame)
end
close(v);
