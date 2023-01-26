function matrixV = formV(Ixi,a)
matrixV = zeros(3,9);

% r1 = Ixi(1);
% r2 = Ixi(2);
% beta_m = Ixi(3);
theta1 = Ixi(4);
theta2 = Ixi(5);

% beta1 = beta_m - theta1;
% beta2 = beta_m + theta2;

% Q = [cos(beta_m) -sin(beta_m);
%     sin(beta_m) cos(beta_m)];
% Q1 = [cos(beta1) -sin(beta1);
%     sin(beta1) cos(beta1)];
% Q2 = [cos(beta2) -sin(beta2);
%     sin(beta2) cos(beta2)];
R1 = [cos(theta1) -sin(theta1);
    sin(theta1) cos(theta1)];
R2 = [cos(theta2) sin(theta2);
    -sin(theta2) cos(theta2)];

% r_G = [-r2 r1];
% r_o1 = [-(r2-a*sin(beta_m)) (r1-a*cos(beta_m))];
% r_o2 = [-(r2+a*sin(beta_m)) (r1+a*cos(beta_m))];

matrixV(1:2,1:6) = [eye(2) R1' R2'];
matrixV(3,1:6) = [0 0 a*sin(theta1) -a*cos(theta1) a*sin(theta2) a*cos(theta2)];
matrixV(3,7:9) = [1 1 1];
end