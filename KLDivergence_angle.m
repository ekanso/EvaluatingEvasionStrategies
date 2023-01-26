function kl = KLDivergence_angle(X,Y,k)
[D,N] = size(X);
[~,M] = size(Y);
rho_p = zeros(N,1);
rho_q = zeros(N,1);
for i = 1:N
    X_relative = [X(:,1:i-1) X(:,i+1:end)] - X(:,i);
    Y_relative = Y - X(:,i);
    if D>1
        X_relative = wrapToPi(X_relative);
        Y_relative = wrapToPi(Y_relative);
        XT_sort = sort(vecnorm(X_relative));
        YT_sort = sort(vecnorm(Y_relative));
    else
        disX = abs(X_relative);
        disY = abs(Y_relative);
%         disX(disX>pi) = 2*pi - disX(disX>pi);
%         disY(disY>pi) = 2*pi - disY(disY>pi);
        XT_sort = sort(disX);
        YT_sort = sort(disY);
    end
    rho_p(i) = XT_sort(k);
    rho_q(i) = YT_sort(k);
    %     rho_p(i) = min(vecnorm(X_trans));
    %     rho_q(i) = min(vecnorm(Y_trans));
end
kl = D/N*sum(-log(rho_p./rho_q)) + log(M/(N-1));
