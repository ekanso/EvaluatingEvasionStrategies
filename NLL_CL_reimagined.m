function out = NLL_CL_reimagined(eta,x)

[D,N] = size(x);

if D ~= 5
    error('Improper data matrix dimension');
end

if length(eta) == 2 || length(eta) == 9
    out = 1/N*sum(-log(pX_CL_reimagined(x,eta)));
else
    error('Invalid parameter vector length');
end

% global indx;

% if isreal(out) && ~isinf(out) && ~isnan(out)
%     
%     if length(eta) == 2
%         fprintf('%i - %f: \t(%f,\t%f)pi\n',indx,out,eta(1)/pi,eta(2)/pi);
%     elseif length(eta) == 9
%         fprintf('%i - %f: \t(%f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f)pi\n',indx,out*N,eta(1)/pi,eta(2)/pi,eta(3)/pi,eta(4)/pi,eta(5)/pi,eta(6)/pi,eta(7)/pi,eta(8)/pi,eta(9)/pi);
%     end
%     
%     figure(99)
%     plot(indx,out,'k.');
%     hold on
%     indx = indx + 1;
%     %ylim([5.2 5.5])
%     
% end

end

