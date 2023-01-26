function matrixM = formM(Ixi,a,a1,a2,m,m1,m2,J,J1,J2)

% beta_m = Ixi(3);
theta1 = Ixi(4);
theta2 = Ixi(5);
% beta1 = beta_m - theta1;
% beta2 = beta_m + theta2;

matrixM = zeros(9,3);

matrixM(1,:) = [(m(1)+m(2)) 0 0];
matrixM(2,:) = [0 (m(1)+m(3)) 0];
matrixM(3,:) = [(m1(1)+m1(2))*cos(theta1) -(m1(1)+m1(2))*sin(theta1) (m1(1)+m1(2))*a*sin(theta1)];
matrixM(4,:) = [(m1(1)+m1(3))*sin(theta1) (m1(1)+m1(3))*cos(theta1) -(m1(1)+m1(3))*(a*cos(theta1)+a1)];
matrixM(5,:) = [(m2(1)+m2(2))*cos(theta2) (m2(1)+m2(2))*sin(theta2) (m2(1)+m2(2))*a*sin(theta2)];
matrixM(6,:) = [-(m2(1)+m2(3))*sin(theta2) (m2(1)+m2(3))*cos(theta2) (m2(1)+m2(3))*(a*cos(theta2)+a2)];
matrixM(7,:) = [0 0 (J(1)+J(2))];
matrixM(8,:) = [-a1*(m1(1)+m1(3))*sin(theta1) -a1*(m1(1)+m1(3))*cos(theta1) (J1(1)+J1(2))+(m1(1)+m1(3))*a*a1*cos(theta1)];
matrixM(9,:) = [-a2*(m2(1)+m2(3))*sin(theta2) a2*(m2(1)+m2(3))*cos(theta2) (J2(1)+J2(2))+(m2(1)+m2(3))*a*a2*cos(theta2)];


end