function p = pX_DO_zero_reimagined(x,eta)

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
    
    p(n)  = p_dp(phi,theta,lam,sig_phi,sig_lam,sig_theta);
end
% toc
end

function p = p_dp(phi,theta,lam,sig_phi,sig_lam,sig_theta)

if (sig_phi<=1e-4) || (sig_theta<=1e-4) || (sig_lam<=1e-4)
    p = 0;
else
%     tic
    a1 = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,pi),0,2*pi,@(phih) -phih,@(phih) pi-phih);
    a2 = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,pi),0,2*pi,@(phih) pi-phih,@(phih) 2*pi-phih);
    p = min(1,max(a1+a2,0));
%     t1=toc;
%     tic
%     a = integral2(@(phih,lamh)intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,pi),0,2*pi,-pi,pi);
%     t2=toc;

end
end


function out = intfunc_alt(phih,lamh,phi,theta,lam,sig_phi,sig_lam,sig_theta,offset)

theta_target = wrapToPi(lamh+phih+offset);
thetah = theta_target - 4/9/pi*abs(theta_target).*theta_target;
thetah(abs(theta_target)>0.98*pi) = (-abs(theta_target(abs(theta_target)>0.98*pi))*254/9/pi+254/9).*theta_target(abs(theta_target)>0.98*pi);
out = VonMises(lamh,lam,sig_lam).*...
    VonMises(phih,phi,sig_phi).*...
    VonMises(thetah,theta,sig_theta);

end
