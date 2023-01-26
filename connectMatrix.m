function A = connectMatrix(alpha1,alpha2,a,m1,m2,J)
Ilock = lockMatrix(alpha1,alpha2,a,m1,m2,J);
Icouple = coupleMatrix(alpha1,alpha2,a,m2,J);

A = -Ilock\Icouple;


function Ilock = lockMatrix(alpha1,alpha2,a,m1,m2,J)
MM = [m1*(1+cos(alpha1)^2+cos(alpha2)^2) + m2*(sin(alpha1)^2+sin(alpha2)^2),...
    (m1-m2)*(sin(2*alpha1)-sin(2*alpha2))/2;
    (m1-m2)*(sin(2*alpha1)-sin(2*alpha2))/2, ...
    m2*(1+cos(alpha1)^2+cos(alpha2)^2) + m1*(sin(alpha1)^2+sin(alpha2)^2)];

HH = [1/2*(m1-m2)*a*(sin(2*alpha1)+sin(2*alpha2)) - m2*a*(sin(alpha1)+sin(alpha2));
    1/2*(m1-m2)*a*(cos(2*alpha2)-cos(2*alpha1)) + m2*a*(cos(alpha1)-cos(alpha2))];

JJ = 3*J + m1*a^2*(sin(alpha1)^2+sin(alpha2)^2)+m2*a^2*((1+cos(alpha1))^2+(1+cos(alpha2))^2);
Ilock = [MM,HH;
    HH',JJ];
end
function Icouple = coupleMatrix(alpha1,alpha2,a,m2,J)
    Icouple = [-m2*a*sin(alpha1), m2*a*sin(alpha2);
        m2*a*cos(alpha1), m2*a*cos(alpha2);
        J+m2*a^2*(1+cos(alpha1)), -J-m2*a^2*(1+cos(alpha2))];
    
end
end