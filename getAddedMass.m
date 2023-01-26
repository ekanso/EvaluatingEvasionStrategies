function [ma1,ma2,ja] = getAddedMass(a,b,c)
factor = 3*a*b*c;
syms lambda
delta = sqrt((a^2+lambda)*(b^2+lambda)*(c^2+lambda));
const_alpha = double(vpaintegral(1/(a^2+lambda)/delta,[0,Inf])*a*b*c);
const_beta = double(vpaintegral(1/(b^2+lambda)/delta,[0,Inf])*a*b*c);

ma1 = const_alpha/(2-const_alpha)*a*b*c/factor;
ma2 = const_beta/(2-const_beta)*a*b*c/factor;
ja = 1/5*(a^2 - b^2)^2*(const_beta - const_alpha)/(2*(a^2 - b^2) + (a^2 + b^2)*(const_alpha - const_beta))*a*b*c/factor;
end