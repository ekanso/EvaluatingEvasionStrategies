function lfit = circular_linearFitting(angle1,angle2)
assert(min(size(angle1)) == 1,'Wrong Input size!!!')
assert(min(size(angle2)) == 1,'Wrong Input size!!!')
assert(isequal(size(angle1), size(angle2)), 'Input sizes do not match!!!')
a = 0.55;
b = -pi/2;
% Gradient descent
% Loss function L = (sin(angle2_fit) - sin(angle2))^2 + (cos(angle2_fit) - cos(angle2))^2
alpha = 0.1; % learning rate
for k = 1:10000
    angle2_fit = a*angle1 + b;
    a_next = a-alpha*mean(sin(angle2_fit - angle2).*angle1);
    b_next = b-alpha*mean(sin(angle2_fit - angle2));
    if abs(a-a_next) + abs(b-b_next) < 1e-6
        break
    end
    a = a_next; b = b_next;
end
if k == 10000
    disp('not converged')
end
lfit = [a,wrapToPi(b)];
end