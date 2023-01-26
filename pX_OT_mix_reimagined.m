function p = pX_OT_mix_reimagined(x,eta)

sig_phi = eta(1);
sig_lam = eta(2);
sig_theta = eta(3);

[D,N] = size(x);

if D~=5
    error('Invalid data size');
end

p = zeros(1,N);

phi_a = x(1,:);
lam_a = x(2,:);

theta_a = x(5,:);
tic
parfor n = 1:N
%     n
    phi = phi_a(n);
    lam = lam_a(n);
    
    theta = theta_a(n);
    
    p(n)  = p_dp(phi,theta,lam,sig_phi,sig_lam,sig_theta);
end
toc
end

function p = p_dp(phi,theta,lam,sig_phi,sig_lam,sig_theta)

if (sig_phi<=1e-4) || (sig_theta<=1e-4) || (sig_lam<=1e-4)
    p = 0;
else
    % 0,pi
    a11 = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,pi/2),0,2*pi,@(phih) max(0,wrapToPi(-phih+pi/2)),@(phih) min(pi,wrapTo2Pi(pi-phih+pi/2)) );    
    a12 = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,pi/2),0,2*pi,@(phih) max(0,wrapToPi(pi-phih+pi/2)),@(phih) min(pi,wrapTo2Pi(2*pi-phih+pi/2)));
    % -pi,0
    a21 = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,3*pi/2),0,2*pi,@(phih) max(pi,wrapTo2Pi(-phih-pi/2))-2*pi,@(phih) min(0,wrapToPi(pi-phih-pi/2)));
    a22 = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,3*pi/2),0,2*pi,@(phih) max(pi,wrapTo2Pi(pi-phih-pi/2))-2*pi,@(phih) min(0,wrapToPi(2*pi-phih-pi/2)));
    p = min(max(a11 + a12 + a21 + a22,0),1);
    
    p = min(1,max(a1_alt + a2_alt,0));
end
end

function out = intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,offset)

theta_target = wrapToPi(lamh+phih+offset);

thetah = (1 - 4/9/pi*abs(theta_target)).*theta_target;
% thetah(abs(theta_target)>0.99*pi) = (-abs(theta_target(abs(theta_target)>0.99*pi))*56/pi+56).*theta_target;
% thetah(abs(theta_target)>0.99*pi) = thetah(abs(theta_target)>0.99*pi) - 
out = VonMises(lamh,lam,sig_lam).*...
    VonMises(phih,phi,sig_phi).*...
    VonMises(thetah,theta,sig_theta);
end
