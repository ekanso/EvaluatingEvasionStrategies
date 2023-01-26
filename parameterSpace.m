clear;
close all;
%% set parameters
a2 = 10;
b2 = 2;
c2 = 10;
a = 10;
b = 2;
c = 10;
a1 = 10;
b1 = 2;
c1 = 10;

L_c = 2*(a + a1 + a2);% to make the size dimensionless
a = a/L_c;
b = b/L_c;
c = c/L_c;
a1 = a1/L_c;
b1 = b1/L_c;
c1 = c1/L_c;
a2 = a2/L_c;
b2 = b2/L_c;
c2 = c2/L_c;

rho = 1000;
M_c = 4/3*rho*(a*b*c+a1*b1*c1+a2*b2*c2);

%% body mass and added mass
m = zeros(3,1);
m1 = zeros(3,1);
m2 = zeros(3,1);
J = zeros(2,1);
J1 = zeros(2,1);
J2 = zeros(2,1);

% calculate the added mass and moment of inertia
[m(2),m(3),m1(2),m1(3),m2(2),m2(3),J(2),J1(2),J2(2)] = getAddedMass(a,b,c,a1,b1,c1,a2,b2,c2);

% change of reference point
J1(2) = J1(2) + a1^2*m1(3);
J2(2) = J2(2) + a2^2*m2(3);
%% compute
rx = 0;ry = 0;
beta_m = 0;
gridNum = 321;
Ax1 = zeros(gridNum);
Ax2 = zeros(gridNum);
Ay1 = zeros(gridNum);
Ay2 = zeros(gridNum);
Ab1 = zeros(gridNum);
Ab2 = zeros(gridNum);
bound = pi;
[x,y] = meshgrid(-bound:bound/(gridNum-1)*2:bound,-bound:bound/(gridNum-1)*2:bound);
for j = 1:gridNum
    for k = 1:gridNum
%         theta_1 = 2*pi*(k-(gridNum+1)/2)/(gridNum-1);theta_2 = 2*pi*(j-(gridNum+1)/2)/(gridNum-1);
        theta_1 = x(j,k);theta_2 = y(j,k);
        Ixi = [rx;ry;beta_m;theta_1;theta_2];
        %% connection method 2
        A = connectMatrix(theta_1,theta_2,a,m(1)+m(2),m(1)+m(3),J(1)+J(2));
        
        Ax1(j,k) = A(1,1);
        Ax2(j,k) = A(1,2);
        Ay1(j,k) = A(2,1);
        Ay2(j,k) = A(2,2);
        Ab1(j,k) = A(3,1);
        Ab2(j,k) = A(3,2);
        %         if (abs(theta_1-0.1)<0.1)&&(abs(theta_2-0.1)<0.1)
        %
        %         end
    end
end
%% plot


fig1 = figure();
fig1.Color = 'w';
xh = quiver(x,y,Ax1,Ax2,'LineWidth',1,'ShowArrowHead','off');
axis equal;
ax1 = gca;title('$\mathbf A_x$ field','Interpreter','latex');
xlabel('$\theta_1$','Interpreter','latex');ylabel('$\theta_2$','Interpreter','latex');
ax1.XTick = -pi:2*pi:pi;
ax1.XTickLabel = {'-$\pi$','$\pi$'};
ax1.YTick = -pi:2*pi:pi;
ax1.YTickLabel = {'-$\pi$','$\pi$'};
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 24; ax1.Box = 'on'; ax1.LineWidth = 2;
xlim([-pi pi]);
ylim([-pi pi]);

fig2 = figure();
fig2.Color = 'w';
yh = quiver(x,y,Ay1,Ay2,'LineWidth',1,'ShowArrowHead','off');
axis equal;
ax2 = gca;title('$\mathbf A_y$ field','Interpreter','latex');
xlabel('$\theta_1$','Interpreter','latex');ylabel('$\theta_2$','Interpreter','latex');
ax2.XTick = -pi:2*pi:pi;
ax2.XTickLabel = {'-$\pi$','$\pi$'};
ax2.YTick = -pi:2*pi:pi;
ax2.YTickLabel = {'-$\pi$','$\pi$'};
ax2.TickLabelInterpreter = 'latex';
ax2.FontSize = 24; ax2.Box = 'on'; ax2.LineWidth = 2;
xlim([-pi pi]);
ylim([-pi pi]);


fig3 = figure();
fig3.Color = 'w';
bh = quiver(x,y,Ab1,Ab2,'LineWidth',1,'ShowArrowHead','off');
axis equal;
ax3 = gca;title('$\mathbf A_\beta$ field','Interpreter','latex');
xlabel('$\theta_1$','Interpreter','latex');ylabel('$\theta_2$','Interpreter','latex');
ax3.XTick = -pi:2*pi:pi;
ax3.XTickLabel = {'-$\pi$','$\pi$'};
ax3.YTick = -pi:2*pi:pi;
ax3.YTickLabel = {'-$\pi$','$\pi$'};
ax3.TickLabelInterpreter = 'latex';
ax3.FontSize = 24; ax3.Box = 'on'; ax3.LineWidth = 2;
xlim([-pi pi]);
ylim([-pi pi]);

% cmap=jet(1e4);
% cmap = gray(1e4);
cmap = zeros(1200,3);
% cmap(1:300,1) = linspace(0.1451,0.1725,300);
% cmap(1:300,2) = linspace(0.2039,0.4980,300);
% cmap(1:300,3) = linspace(0.5804,0.7216,300);
% cmap(300:600,1) = linspace(0.1725,0.2549,301);
% cmap(300:600,2) = linspace(0.4980,0.7137,301);
% cmap(300:600,3) = linspace(0.7216,0.7686,301);
% cmap(600:900,1) = linspace(0.2549,0.6314,301);
% cmap(600:900,2) = linspace(0.7137,0.8549,301);
% cmap(600:900,3) = linspace(0.7686,0.7059,301);
% cmap(900:1200,1) = linspace(0.6314,1,301);
% cmap(900:1200,2) = linspace(0.8549,1,301);
% cmap(900:1200,3) = linspace(0.7059,0.8,301);

% red, white, blue
cmap(1:600,1) = linspace(0,1,600).^1;
cmap(1:600,2) = linspace(0.6,1,600).^1.5;
cmap(1:600,3) = linspace(0.9,1,600).^1.5;

cmap(600:1200,1) = linspace(1,0.9,601).^1.5;
cmap(600:1200,2) = linspace(1,0.6,601).^1.5;
cmap(600:1200,3) = linspace(1,0,601).^1;

[curlx,~]= curl(x,y,Ax1,Ax2);
% [curlx2,~]= curl(x,y,Ax21,Ax22);
curlx_max = max((-curlx(:)));
curlx = curlx/curlx_max;
fig4 = figure;
fig4.Color = 'w';
pcolor(x,y,curlx);
shading(gca,'interp');
% contourf(x,y,curlz,'ShowText','off','EdgeColor','none');
view(2);
colormap(cmap);
colorbar('Ticks',[-1,0,1]);
caxis([-1 1]);
axis equal;
ax4 = gca;
% title(['Curl($\mathbf A_x$), Scale=' num2str(curlx_max,3)],'Interpreter','latex');
title('Curl($\mathbf A_x$)','Interpreter','latex');
xlabel('$\alpha_1$','Interpreter','latex');ylabel('$\alpha_2$','Interpreter','latex');
ax4.XTick = -pi:2*pi:pi;
ax4.XTickLabel = {'-$\pi$','$\pi$'};
ax4.YTick = -pi:2*pi:pi;
ax4.YTickLabel = {'-$\pi$','$\pi$'};
ax4.TickLabelInterpreter = 'latex';
ax4.FontSize = 24; ax4.Box = 'on'; ax4.LineWidth = 2;
xlim([-pi pi]);
ylim([-pi pi]);hold on;

[curly,~]= curl(x,y,Ay1,Ay2);
curly_max = max(abs(curly(:)));
curly = curly/curly_max;
fig5 = figure;
fig5.Color = 'w';
pcolor(x,y,curly);
shading interp;
% contourf(x,y,curlz,'ShowText','off','EdgeColor','none');
view(2);
colormap(cmap);
colorbar('Ticks',[-1,0,1]);
caxis([-1 1]);
axis equal;
ax5 = gca;
% title(['Curl($\mathbf A_y$), Scale=' num2str(curly_max,3)],'Interpreter','latex');
title('Curl($\mathbf A_y$)','Interpreter','latex');
xlabel('$\alpha_1$','Interpreter','latex');ylabel('$\alpha_2$','Interpreter','latex');
ax5.XTick = -pi:2*pi:pi;
ax5.XTickLabel = {'-$\pi$','$\pi$'};
ax5.YTick = -pi:2*pi:pi;
ax5.YTickLabel = {'-$\pi$','$\pi$'};
ax5.TickLabelInterpreter = 'latex';
ax5.FontSize = 24; ax5.Box = 'on'; ax5.LineWidth = 2;
xlim([-pi pi]);
ylim([-pi pi]);hold on;


[curlb,~]= curl(x,y,Ab1,Ab2);
curlb_max = max(abs(curlb(:)));
curlb = curlb/curlb_max;
fig6 = figure;
fig6.Color = 'w';
pcolor(x,y,curlb);
shading interp;
% contourf(x,y,curlz,'ShowText','off','EdgeColor','none');
view(2);
colormap(cmap);
colorbar('Ticks',[-1,0,1]);
caxis([-1 1]);
axis equal;
ax6 = gca;
% title(['Curl($\mathbf A_\beta$), Scale=' num2str(curlb_max,3)],'Interpreter','latex');
title('Curl($\mathbf A_\beta$)','Interpreter','latex');
xlabel('$\alpha_1$','Interpreter','latex');ylabel('$\alpha_2$','Interpreter','latex');
ax6.XTick = -pi:2*pi:pi;
ax6.XTickLabel = {'-$\pi$','$\pi$'};
ax6.YTick = -pi:2*pi:pi;
ax6.YTickLabel = {'-$\pi$','$\pi$'};
ax6.TickLabelInterpreter = 'latex';
ax6.FontSize = 24; ax6.Box = 'on'; ax6.LineWidth = 2;
xlim([-pi pi]);
ylim([-pi pi]);
%% generate path
hold on;
load('data/exp_data_Cstart.mat');
theta1 = theta1/180*pi;
theta2 = theta2/180*pi;
plot(theta1,theta2,'.r');

time = 0:pi/100:2*pi;
a0 =       17.94;
a1 =      -1.059;
b1 =       23.23;
a2 =       -29.4;
b2 =       30.49;
a3 =       14.93;
b3 =      -31.43;
w =      0.7886;
theta2 = a0 + a1*cos(time*w) + b1*sin(time*w) + a2*cos(2*time*w) + b2*sin(2*time*w) + a3*cos(3*time*w) + b3*sin(3*time*w);
theta2 = theta2/180*pi;
a0 =       13.46;
a1 =       2.326;
b1 =       8.468;
a2 =       -50.6;
b2 =      -27.37;
a3 =       33.43;
b3 =       23.41;
w =      0.7615;
theta1 = a0 + a1*cos(time*w) + b1*sin(time*w) + a2*cos(2*time*w) + b2*sin(2*time*w) + a3*cos(3*time*w) + b3*sin(3*time*w);
theta1 = theta1/180*pi;


plot(ax6,theta1,theta2,'k','LineWidth',2);
