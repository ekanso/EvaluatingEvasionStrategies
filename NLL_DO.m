function out = NLL_DO(eta,x,chi)

[D,N] = size(x);

if D ~= 5
    error('Improper data matrix dimension');
end

% out = 0;
if (length(eta) == 2) && (chi~=0)
    eta(3) = eta(1);
    out = 1/N*sum(-log(pX_DO_mix(x,eta,chi)));
elseif (length(eta) == 1) && (chi==0)
    eta(2) = eta(1);
    eta(3) = eta(1);
    out = 1/N*sum(-log(pX_DO_zero(x,eta)));
% elseif length(eta) == 4
%     out = 1/N*sum(-log(pX_S4_mix(x,eta(1:3),eta(4))))
else
    error('Invalid parameter vector length');
end


end

