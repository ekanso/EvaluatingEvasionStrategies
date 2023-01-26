function matrixM = formM_old(Ixi,a,a1,a2,m,m1,m2,J,J1,J2)

beta_m = Ixi(3);
theta1 = Ixi(4);
theta2 = Ixi(5);
beta1 = beta_m - theta1;
beta2 = beta_m + theta2;

matrixM = zeros(9,5);

matrixM(1,:) = [(m(1)+m(2))*cos(beta_m) (m(1)+m(2))*sin(beta_m) 0 0 0];
matrixM(2,:) = [-(m(1)+m(3))*sin(beta_m) (m(1)+m(3))*cos(beta_m) 0 0 0];
matrixM(3,:) = [(m1(1)+m1(2))*cos(beta1) (m1(1)+m1(2))*sin(beta1) (m1(1)+m1(2))*a*sin(theta1) 0 0];
matrixM(4,:) = [-(m1(1)+m1(3))*sin(beta1) (m1(1)+m1(3))*cos(beta1) -(m1(1)+m1(3))*(a*cos(theta1)+a1) (m1(1)+m1(3))*a1 0];
matrixM(5,:) = [(m2(1)+m2(2))*cos(beta2) (m2(1)+m2(2))*sin(beta2) (m2(1)+m2(2))*a*sin(theta2) 0 0];
matrixM(6,:) = [-(m2(1)+m2(3))*sin(beta2) (m2(1)+m2(3))*cos(beta2) (m2(1)+m2(3))*(a*cos(theta2)+a2) 0 (m2(1)+m2(3))*a2];
matrixM(7,:) = [0 0 (J(1)+J(2)) 0 0];
matrixM(8,:) = [a1*(m1(1)+m1(3))*sin(beta1) -a1*(m1(1)+m1(3))*cos(beta1) (J1(1)+J1(2))+(m1(1)+m1(3))*a*a1*cos(theta1) -(J1(1)+J1(2)) 0];
matrixM(9,:) = [-a2*(m2(1)+m2(3))*sin(beta2) a2*(m2(1)+m2(3))*cos(beta2) (J2(1)+J2(2))+(m2(1)+m2(3))*a*a2*cos(theta2) 0 (J2(1)+J2(2))];


end