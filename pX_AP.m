function p = pX_AP(x,eta)

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
if (sig_phi<=1e-4) || (sig_theta<=1e-4)
    p = 0;
% elseif (sig_phi>=0.08) && (sig_theta>=0.08)
%     p = integral(@(thetah)integrand(thetah,theta,phi,sig_phi,sig_theta),0,2*pi);
% elseif (sig_phi>=0.08)
%     p = integral(@(thetah)integrand4(thetah,theta,phi,sig_phi,sig_theta),0,2*pi);
% elseif (sig_theta>=0.08)
%     p = integral(@(thetah)integrand5(thetah,theta,phi,sig_phi,sig_theta),0,2*pi);
else
    p = integral(@(thetah)integrand3(thetah,theta,phi,sig_phi,sig_theta),0,2*pi);
end
end

function out = integrand(thetah,theta,phi,sig_phi,sig_theta)

arg = cos(thetah + pi - phi)/sig_phi^2+cos(thetah-theta)/sig_theta^2;
out = exp(arg)/(4*pi^2*besseli(0,1/sig_phi^2)*besseli(0,1/sig_theta^2));

end

function out = integrand2(thetah,theta,phi,sig_phi,sig_theta)
% first order approximation
arg = (cos(thetah + pi - phi)-1)/sig_phi^2 + (cos(thetah-theta)-1)/sig_theta^2;
out = exp(arg)/(2*pi*sig_theta*sig_phi);

end

function out = integrand3(thetah,theta,phi,sig_phi,sig_theta)
% second order approximation

out = VonMises(thetah+pi,phi,sig_phi).*VonMises(thetah,theta,sig_theta);

end

function out = integrand4(thetah,theta,phi,sig_phi,sig_theta)
% second order approximation for only sig_theta
arg = cos(thetah + pi - phi)/sig_phi^2 + (cos(thetah-theta)-1)/sig_theta^2;
out = exp(arg)/(4*pi^2*besseli(0,1/sig_phi^2)*(sig_theta+sig_theta^3/8)/sqrt(2*pi));

end

function out = integrand5(thetah,theta,phi,sig_phi,sig_theta)
% second order approximation for only sig_phi
arg = (cos(thetah + pi - phi)-1)/sig_phi^2 + cos(thetah-theta)/sig_theta^2;
out = exp(arg)/(4*pi^2*(sig_phi+sig_phi^3/8)/sqrt(2*pi)*besseli(0,1/sig_theta^2));

end
