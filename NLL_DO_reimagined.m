function out = NLL_DO_reimagined(eta,x,chi)

[D,N] = size(x);

if D ~= 5
    error('Improper data matrix dimension');
end

% out = 0;
if length(eta) == 3
    if chi == 0
    out = 1/N*sum(-log(pX_DO_zero_reimagined(x,eta)));    
    else
    out = 1/N*sum(-log(pX_DO_mix_reimagined(x,eta,chi)));
    end
% elseif length(eta) == 4
%     out = 1/N*sum(-log(pX_S4_mix(x,eta(1:3),eta(4))))
else
    error('Invalid parameter vector length');
end

% disp(['DO: ','chi=',num2str(chi),'eta=',num2str(eta'),'NLL=',num2str(out)])
end

