function out = NLL_CL(eta,x)

[D,N] = size(x);

if D ~= 5
    error('Improper data matrix dimension');
end

if length(eta) == 2 || length(eta) == 9
    out = 1/N*sum(-log(pX_CL(x,eta)));
else
    error('Invalid parameter vector length');
end

end

