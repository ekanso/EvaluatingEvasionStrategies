function p = pX_DO_mix_reimagined(x,eta,chi)

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
% eta
% tic
% ticBytes(gcp)
parfor n = 1:N
%     n
    phi = phi_a(n);
    lam = lam_a(n);
    theta = theta_a(n);
    
    p(n)  = p_dp(phi,theta,lam,sig_phi,sig_lam,sig_theta,chi);
end
% tocBytes(gcp)
% toc
end
% 1.0939,0.6608,0.1663
function p = p_dp(phi,theta,lam,sig_phi,sig_lam,sig_theta,chi)

if (sig_phi<=1e-4) || (sig_theta<=1e-4) || (sig_lam<=1e-4)
    p = 0;
else
    % 0,pi
%     tic
    a111 = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,pi-chi),chi,pi+chi, 0 ,@(phih) pi-phih+chi, 'RelTol',1e-6);
    a112 = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,pi-chi),pi+chi,2*pi+chi,@(phih) 2*pi-phih+chi,pi, 'RelTol',1e-6);
    
    a121 = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,pi-chi),chi,pi+chi,@(phih) pi-phih+chi, pi, 'RelTol',1e-6);
    a122 = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,pi-chi),pi+chi,2*pi+chi,0,@(phih) 2*pi-phih+chi, 'RelTol',1e-6);
%     t1 = toc;
%     tic
%     a1_alt = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,pi-chi),0,2*pi,0,pi);
%     t2 = toc;
    % -pi,0
%     tic
    a211 = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,pi+chi),-chi,pi-chi,@(phih) -phih-chi,0, 'RelTol',1e-6);
    a212 = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,pi+chi),pi-chi,2*pi-chi,-pi,@(phih) pi-phih-chi, 'RelTol',1e-6);

    a221 = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,pi+chi),-chi,pi-chi, -pi,@(phih) -phih-chi, 'RelTol',1e-6);
    a222 = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,pi+chi),pi-chi,2*pi-chi,@(phih) pi-phih-chi,0, 'RelTol',1e-6);
%     t3 = toc;
%     tic
%     a2_alt = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,pi+chi),0,2*pi,-pi,0);
%     t4 = toc;
    p = min(max(a111 + a112 + a121 + a122 + a211 + a212 + a221 + a222,0),1);
%     p2 = a1_alt+a2_alt
end
end

function out = intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,offset)

theta_target = wrapToPi(lamh+phih+offset);
thetah = theta_target - 4/9/pi*abs(theta_target).*theta_target;
out = VonMises(lamh,lam,sig_lam).*...
    VonMises(phih,phi,sig_phi).*...
    VonMises(thetah,theta,sig_theta);
end
