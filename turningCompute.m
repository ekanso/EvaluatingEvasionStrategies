function turning = turningCompute(aa,bb)

frame_num = 180;
time_expect = 2*pi;% total simulation time
mode = 3;

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

%% body mass and added mass
m = zeros(3,1);
J = zeros(2,1);
J1 = zeros(2,1);
J2 = zeros(2,1);
m(1) = a*b*c/(a*b*c + a1*b1*c1 + a2*b2*c2);
J(1) = 1/5*(a^2+b^2)*m(1);
J1(1) = (1/5*(a1^2+b1^2)+a1^2)*m(1);
J2(1) = (1/5*(a2^2+b2^2)+a2^2)*m(1);
% calculate the added mass and moment of inertia
[m(2),m(3),J(2)] = getAddedMass(a,b,c);

% change of reference point
J1(2) = J1(2) + a1^2*m(3);
J2(2) = J2(2) + a2^2*m(3);

%% initialization of the variations
[xi,Ixi] = initialConfig(mode,1,aa,bb);
% total momentums

beta_m = Ixi(3);
theta1 = Ixi(4);
theta2 = Ixi(5);

beta_mean_ini = beta_m - theta1/3 + theta2/3;

%% time integration
var = Ixi(1:3);

Momt = [0;0;0];

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,var] = ode45(@(t,var) firstorderdt(t,var,Momt,a,a1,a2,m,m,m,J,J1,J2,aa,bb),linspace(0,time_expect,frame_num),var,opts);
[theta1,theta2,Dtheta1,Dtheta2] = prescribedAngle(t(end),0,aa,bb);

beta_m = var(end,3);

beta_mean = beta_m-theta1/3+theta2/3-beta_mean_ini;

turning  = -beta_mean;

