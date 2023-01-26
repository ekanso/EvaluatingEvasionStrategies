function p = pX_CL_reimagined(x,eta)

sig_phi = eta(1);
sig_theta = eta(2);


[D,N] = size(x);

if D~=5
    error('Invalid data size');
end

p = zeros(1,N);

phi_a = x(1,:);

theta_a = x(5,:);

parfor n = 1:N
    phi = phi_a(n);
    theta = theta_a(n);
    
    p(n)  = p_dp(phi,theta,sig_phi,sig_theta);
end
end



function p = p_dp(phi,theta,sig_phi,sig_theta)

% phid = xi(phi,sig_phi);
phid = min(integral(@(phih)VonMises(phih,phi,sig_phi),0,pi),1);
% if (sig_theta>= 0.08)
%     p = 1/(2*pi*besseli(0,1/sig_theta^2))*(exp(cos(theta+pi/2)/sig_theta^2)*phid+...
%         exp(cos(theta-pi/2)/sig_theta^2)*(1-phid));
% elseif (sig_theta >= 1e-8)
%     p = 1/(sqrt(2*pi)*(sig_theta+sig_theta^3/8))*(exp((cos(theta+pi/2)-1)/sig_theta^2)*phid+...
%         exp((cos(theta-pi/2)-1)/sig_theta^2)*(1-phid));
% end

p = phid*VonMises(theta,-pi/18*7,sig_theta) + (1-phid)*VonMises(theta,pi/18*7,sig_theta);
end

function out = xi(phi,sig_phi)

if (sig_phi >= 1e-5)
    out = integral(@(phih)VonMises(phih,phi,sig_phi),0,pi);
else
    out = ((phi>0)&&(phi<pi))*1 + (phi==0)*0.5 + (phi==pi)*0.5;
end
out(out>1) = 1;
end
