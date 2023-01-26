clear;
close all;
set(groot,'defaultLineLineWidth',2);
set(groot,'defaultFigureColor','w');
set(groot,'defaultTextFontsize',22);
set(groot,'defaultAxesFontsize',22);
set(groot,'defaultPolarAxesFontsize',22);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultPolarAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesLineWidth',1);
%%
% bb is major axis, aa is minor axis
aa = zeros(numel(0.05:0.05:2.1),1);
turn = zeros(numel(0.05:0.05:2.1),1);
k = 0;

for bb = 0.05:0.05:2.1
    k = k+1;
    a0 = 0.2+bb/4;
    
    A = -1;
    
    options = optimoptions('fmincon','Display','notify-detailed');
    [aa(k),turn(k)] = fmincon(@(aa)turningCompute(aa,bb),a0,[],[],[],[],0.01,max(2,bb),[],options);
end
bb = 0.05:0.05:2.1;
F = griddedInterpolant(bb*sqrt(2),abs(turn));
%%
% load('result_wmass.mat');
load('result.mat');
xx = zeros(14,121);
yy = zeros(14,121);
bb = majoraxis;
aa = minoraxis;
for k = 5:8:42
    b = 0.05*k;
    t = linspace(0,2*pi,121);
    xx(k,:) = b*(cos(t)+1);
    yy(k,:) = aa(k)*sin(t);
end
example = (5:8:42)*0.05*sqrt(2);

% open('betaMap_wmass.fig');hold on;
open('betaMap.fig');hold on;
% plot(xx',yy','k');
plot(xx'*cos(pi/4)-yy'*sin(pi/4),yy'*cos(pi/4)+xx'*sin(pi/4),'k','linewidth',1.5);
n=0;
for e = example
plot(e,e,'*','color','r');
n = n + 1;
end
turn = F(majoraxis*sqrt(2));

figure,hold on;
plot(majoraxis*sqrt(2),turn);

plot(example,F(example),'r*');
xlabel('max curling angle $\alpha_{\rm{max}}$');
ylabel('turn angle $\theta$');
xlim([0,pi]);
xticks([0,pi/2,pi]);xticklabels({'0','$\pi/2$','$\pi$'});
yticks([0,pi/2,pi]);yticklabels({'0','$\pi/2$','$\pi$'});
% axis equal;
%%
a2 = 1;
b2 = 0.2;
c2 = 1;
a = 1;
b = 0.2;
c = 1;
a1 = 1;
b1 = 0.2;
c1 = 1;

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
M_c = 4/3*rho*(a*b*c+a1*b1*c1+a2*b2*c2)*pi;

m = zeros(3,1);
m1 = zeros(3,1);
m2 = zeros(3,1);
J = zeros(2,1);
J1 = zeros(2,1);
J2 = zeros(2,1);
% m(1) = a*b*c/(a*b*c + a1*b1*c1 + a2*b2*c2);
% m1(1) = a1*b1*c1/(a*b*c + a1*b1*c1 + a2*b2*c2);
% m2(1) = a2*b2*c2/(a*b*c + a1*b1*c1 + a2*b2*c2);
% J(1) = 1/5*(a^2+b^2)*m(1);
% J1(1) = (1/5*(a1^2+b1^2)+a1^2)*m1(1);
% J2(1) = (1/5*(a2^2+b2^2)+a2^2)*m2(1);
% calculate the added mass and moment of inertia
[m(2),m(3),m1(2),m1(3),m2(2),m2(3),J(2),J1(2),J2(2)] = getAddedMass(a,b,c,a1,b1,c1,a2,b2,c2);

% change of reference point
J1(2) = J1(2) + a1^2*m1(3);
J2(2) = J2(2) + a2^2*m2(3);

load('result.mat');
% load('result_wmass.mat');
turn = F(majoraxis*sqrt(2));
total_energyRL = zeros(42,1);
total_energySimple = zeros(42,1);
for k = 1:42
    b = majoraxis(k);
    a = minoraxis(k);
    % fun_theta1 = @(t) abs(b*sin(t));
    % fun_theta2 = @(t) abs(a*cos(t));
    %
    % fun_s = @(t) sqrt(b^2*sin(t)^2 + a^2*cos(t)^2);
    %
    %
    % t0 = 0;
    % t1 = fzero(@(te) (integral(fun_theta1,t0,te)-0.1), t0);
    % t2 = fzero(@(te) (integral(fun_theta2,t0,te)-0.1), t0);
    
    % fun_energy = (dalpha1^2+dalpha2^2)*J1(2);
    total_energyRL(k) = J1(2)*integral(@(t) fun_energy_RL(a,b,t),0,2*pi);
    total_energySimple(k) = J1(2)*integral(@(t) fun_energy_simple(a,b,t),0,2*pi);
end
figure,hold on;
plot(majoraxis*sqrt(2),total_energyRL);
plot(majoraxis*sqrt(2),total_energySimple);
xlabel('max curling angle $\alpha$');
ylabel('effort');
xlim([0,pi]);
xticks([0,pi/2,pi]);xticklabels({'0','$\pi/2$','$\pi$'});
legend('limit on angular velocity','fixed time','location','northwest');

figure,hold on;
plot(abs(turn),total_energyRL);
plot(abs(turn),total_energySimple);
xlabel('turn angle $\theta$');
ylabel('effort');
xlim([0,pi]);
xticks([0,pi/2,pi]);xticklabels({'0','$\pi/2$','$\pi$'});
legend('limit on angular velocity','fixed time','location','northwest');

%%
T = readtable('voesenek2019headtailangleVSturn.csv');
headtotailAng = T.Var1;
TurnAng = T.Var2;
maxcurling = headtotailAng/2;
figure();
plot(maxcurling,TurnAng,'ko');hold on;
load('result.mat');

maxcurling = sort(maxcurling);
turnFromModel = F(maxcurling/180*pi)*180/pi;
plot(maxcurling,turnFromModel,'*-','color','k');
load('result_wmass.mat');
turnFromModel = F(maxcurling/180*pi)*180/pi;
plot(maxcurling,turnFromModel,'<-.','color','k','markerfacecolor','k','markeredgecolor','none','markersize',8);
xlabel('maximum curling angle [deg]');
ylabel('turn angle [deg]');
xlim([25,125]);
ylim([0,120]);
% fitobject = fit(maxcurling,TurnAng,'power2');
% xaxis = 25/180*pi:pi/180:125/180*pi;
% plot(xaxis,fitobject.a*xaxis.^fitobject.b+fitobject.c,'k');
xticks(25:50:125); xticklabels({'50','150','250'});
yticks(0:60:120); yticklabels({'0','60','120'});
legend('Voesenek2019','massless','neutually buoyant','location','northwest');
%%
function energy = fun_energy_RL(a,b,t)
dtdtau = 1./max(abs(-b*sin(t)),abs(a*cos(t)));
energy = (1 + (dtdtau.*min(abs(-b*sin(t)),abs(a*cos(t)))).^2)./dtdtau;
end

function energy = fun_energy_simple(a,b,t)
energy = (-b*sin(t)).^2+(a*cos(t)).^2;
end