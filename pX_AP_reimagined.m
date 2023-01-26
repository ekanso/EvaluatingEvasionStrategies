function p = pX_AP_reimagined(x,eta)

sig_phi = eta(1);
sig_theta = eta(2);


[D,N] = size(x);

if D~=5
    error('Invalid data size');
end

p = zeros(1,N);

phi_a = x(1,:);

theta_a = x(5,:);

for n = 1:N
    phi = phi_a(n);
    theta = theta_a(n);
    
    p(n)  = p_dp(phi,theta,sig_phi,sig_theta);
end
end




function p = p_dp(phi,theta,sig_phi,sig_theta)
if (sig_phi<=1e-4) || (sig_theta<=1e-4)
    p = 0;
else
    p = integral(@(phih)integrand(phih,theta,phi,sig_phi,sig_theta),0,2*pi);
    p = min(max(p,0),1);
end
end


function out = integrand(phih,theta,phi,sig_phi,sig_theta)
% second order approximation
theta_target = wrapToPi(phih - pi);
thetah = theta_target - 4/9/pi*abs(theta_target).*theta_target;
out = VonMises(phih,phi,sig_phi).*VonMises(theta,thetah,sig_theta);
end
