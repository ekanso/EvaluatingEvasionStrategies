function vel = threelink_dynamics(t,act_mode,var,a,m,J,ii)

[alpha1,alpha2,Dalpha1,Dalpha2] = prescribedAngle(t,act_mode,ii);

beta_m = var(3);
Q = [cos(beta_m) -sin(beta_m);
    sin(beta_m) cos(beta_m)];

A = connectMatrix(alpha1,alpha2,a,m(1)+m(2),m(1)+m(3),J(1)+J(2));
vel = [Q zeros(2,1);0 0 1]*A*[Dalpha1;Dalpha2];
end
