function [xi_int,Ixi_int] = initialConfig(act_mode,ii)
xi_int = zeros(5,1);
Ixi_int = zeros(5,1);

%% theta1,theta2 both prescribed
[Ixi_int(4),Ixi_int(5),xi_int(4),xi_int(5)] = prescribedAngle(0,act_mode,ii);

end