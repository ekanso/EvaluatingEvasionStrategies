function out = NLL_AP(eta,x)

[D,N] = size(x);

if D ~= 5
    error('Improper data matrix dimension');
end

if length(eta) == 1
    eta(2) = eta(1);
    out = 1/N*sum(-log(pX_AP(x,eta)));
else
    error('Invalid parameter vector length');
end


end

