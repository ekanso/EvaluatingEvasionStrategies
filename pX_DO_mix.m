function p = pX_DO_mix(x,eta,chi)

sig_phi = eta(1);
sig_lam = eta(2);
sig_theta = eta(3);

[D,N] = size(x);

if D~=5
    error('Invalid data size');
end

p = zeros(1,N);
%     p_phi   = 1/(2*pi);
%     p_lam   = 1/(2*pi);
%     p_u     = 1/(0.20);
%     p_v     = 1/(0.20);
phi_a = x(1,:);
lam_a = x(2,:);
theta_a = x(5,:);
% tic
parfor n = 1:N
%     n
    phi = phi_a(n);
    lam = lam_a(n);
    
    theta = theta_a(n);
    
    p(n)  = p_dp(phi,theta,lam,sig_phi,sig_lam,sig_theta,chi);
end
% toc
end

function p = p_dp(phi,theta,lam,sig_phi,sig_lam,sig_theta,chi)

if (sig_phi<=1e-4) || (sig_theta<=1e-4) || (sig_lam<=1e-4)
    p = 0;
else
    a1 = integral2(@(thetah,lamh)intfunc(thetah,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,-pi+chi),0,2*pi,0,pi);
    a2 = integral2(@(thetah,lamh)intfunc(thetah,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,-pi-chi),0,2*pi,-pi,0);
    a1_alt = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,pi-chi),0,2*pi,0,pi);
    a2_alt = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,pi+chi),0,2*pi,-pi,0);
    
    p = min(1,max((a1 + a2 + a1_alt + a2_alt)/2,0));
end
end

function out = intfunc(thetah,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,offset)

out = VonMises(lamh,lam,sig_lam).*...
    VonMises(thetah-lamh+offset,phi,sig_phi).*...
    VonMises(thetah,theta,sig_theta);

end
function out = intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,offset)

out = VonMises(lamh,lam,sig_lam).*...
    VonMises(phih,phi,sig_phi).*...
    VonMises(lamh+phih+offset,theta,sig_theta);

end