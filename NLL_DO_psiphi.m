function out = NLL_DO_psiphi(eta,x,chi)

[D,N] = size(x);

if D ~= 5
    error('Improper data matrix dimension');
end

% out = 0;
if (length(eta) == 3) && (chi~=0)
    out = 1/N*sum(-log(pX_DO_psiphi(x,eta,chi)));
else
    error('Invalid parameter vector length');
end


end

