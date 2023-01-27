function [alpha1,alpha2,Dalpha1,Dalpha2] = prescribedAngle(time,act_mode,aa,bb)
% theta1 = 0;
% theta2 = 0;
% Dtheta1 = 0;
% Dtheta2 = 0;
r = pi/3;
switch act_mode
    case 'SimpleSwim'
        phase = time + pi/4;
        alpha1 =  cos(phase)*r;
        alpha2 =  sin(phase)*r;
        Dalpha1 = -sin(phase)*r;
        Dalpha2 = cos(phase)*r;
    case 'SimpleTurn'
        phase = time - 3*pi/4;
        alpha1 =  (1/sqrt(2)+cos(phase))*r;
        alpha2 =  (1/sqrt(2)+sin(phase))*r;
        Dalpha1 = (-sin(phase))*r;
        Dalpha2 =  (cos(phase))*r;
    case 'TurnPlusSwim'
        if time<=5*pi/3
            phase = time - 3*pi/4;
            alpha1 =  (1/sqrt(2)+cos(phase))*r;
            alpha2 =  (1/sqrt(2)+sin(phase))*r;
            Dalpha1 = (-sin(phase))*r;
            Dalpha2 =  (cos(phase))*r;
        else
            phase = time - 5*pi/3 + 7*pi/12;
            alpha1 =  cos(phase)*r;
            alpha2 =  sin(phase)*r;
            Dalpha1 = -sin(phase)*r;
            Dalpha2 = cos(phase)*r;
        end
    case 'CStartFit'
        a0 =       17.94;
        a1 =      -1.059;
        b1 =       23.23;
        a2 =       -29.4;
        b2 =       30.49;
        a3 =       14.93;
        b3 =      -31.43;
        w =      0.7886;

        alpha1 = a0 + a1*cos(time*w) + b1*sin(time*w) + a2*cos(2*time*w) + b2*sin(2*time*w) + a3*cos(3*time*w) + b3*sin(3*time*w);
        Dalpha1 = -w*a1*sin(time*w) + w*b1*cos(time*w) - 2*w*a2*sin(2*time*w) + 2*w*b2*cos(2*time*w) - 3*w*a3*sin(3*time*w) + 3*w*b3*cos(3*time*w);
        alpha1 = alpha1/180*pi;
        Dalpha1 = Dalpha1/180*pi;

        a0 =       13.46;
        a1 =       2.326;
        b1 =       8.468;
        a2 =       -50.6;
        b2 =      -27.37;
        a3 =       33.43;
        b3 =       23.41;
        w =      0.7615;

        alpha2 = a0 + a1*cos(time*w) + b1*sin(time*w) + a2*cos(2*time*w) + b2*sin(2*time*w) + a3*cos(3*time*w) + b3*sin(3*time*w);
        Dalpha2 = -w*a1*sin(time*w) + w*b1*cos(time*w) - 2*w*a2*sin(2*time*w) + 2*w*b2*cos(2*time*w) - 3*w*a3*sin(3*time*w) + 3*w*b3*cos(3*time*w);
        alpha2 = alpha2/180*pi;
        Dalpha2 = Dalpha2/180*pi;
    case 'CStartFamily'
        assert(nargin == 4,'argument of axis length needed')

        phi = 3*pi/4;

            alpha1 =  cos(phi)*cos(time+pi/2)*aa - sin(phi)*sin(time+pi/2)*bb + sqrt(2)*bb/2;
            alpha2 =  sin(phi)*cos(time+pi/2)*aa + cos(phi)*sin(time+pi/2)*bb + sqrt(2)*bb/2;
            Dalpha1 = -cos(phi)*sin(time+pi/2)*aa - sin(phi)*cos(time+pi/2)*bb;
            Dalpha2 = -sin(phi)*sin(time+pi/2)*aa + cos(phi)*cos(time+pi/2)*bb;
        
    otherwise
        warning('unknown activation mode')
end


%% ellipse about x = y
% ll = bb;
% bb = 1.2;


%%
% time = mod(time,5/8*pi);
% if time<=pi/2
%     theta1 = -(1-cos(2*time))*cos(time)*pi;
%     theta2 = -(1-cos(2*time))*sin(time);
%     Dtheta1 = -2*sin(2*time)*cos(time)*pi + (1-cos(2*time))*sin(time)*pi;
%     Dtheta2 = -2*sin(2*time)*sin(time) - (1-cos(2*time))*cos(time);
% elseif time<=5/8*pi
%     theta1 = -(1-cos(8*time+3*pi))*cos(time);
%     theta2 = -(1-cos(8*time+3*pi))*sin(time);
%     Dtheta1 = -8*sin(8*time+3*pi)*cos(time) + (1-cos(8*time+3*pi))*sin(time);
%     Dtheta2 = -8*sin(8*time+3*pi)*sin(time) - (1-cos(8*time+3*pi))*cos(time);
% else
%     theta1 = 0;
%     theta2 = 0;
%     Dtheta1 = 0;
%     Dtheta2 = 0;
% end


end
