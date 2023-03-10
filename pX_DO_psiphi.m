function p = pX_DO_psiphi(x,eta,chi)

sig_phi = eta(1);
sig_psi = eta(2);
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
psi_a = x(2,:);
theta_a = x(5,:);
% tic
parfor n = 1:N
%     n
    phi = phi_a(n);
    psi = psi_a(n);
    
    theta = theta_a(n);
    
    p(n)  = p_dp(phi,theta,psi,sig_phi,sig_psi,sig_theta,chi);
end
% toc
end

function p = p_dp(phi,theta,psi,sig_phi,sig_psi,sig_theta,chi)

if (sig_phi<=1e-4) || (sig_theta<=1e-4) || (sig_psi<=1e-4)
    p = 0;
else
    intfunc1 = @(thetah, phih) VonMises(phih,phi,sig_phi).*...
        VonMises(thetah-chi,psi,sig_psi).*...
        VonMises(thetah,theta,sig_theta);
    intfunc2 = @(thetah, phih) VonMises(phih,phi,sig_phi).*...
        VonMises(thetah+chi,psi,sig_psi).*...
        VonMises(thetah,theta,sig_theta);

    a1 = integral2(intfunc1,0,2*pi, @(thetah) thetah-chi-pi, @(thetah) thetah-chi);
    a2 = integral2(intfunc2,0,2*pi, @(thetah) thetah+chi, @(thetah) thetah+chi+pi);
    p = min(1,max(a1+a2,0));
end
end