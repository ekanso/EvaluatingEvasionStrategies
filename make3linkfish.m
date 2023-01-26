function h_fish = make3linkfish(a,b,c,rx,ry,beta_m,alpha1,alpha2,draw_mode,h_fish)
% disp(nargin)
beta1 = beta_m+alpha1;
beta2 = beta_m-alpha2;
r2x = rx - cos(beta_m)*a - cos(beta_m-alpha2)*a;
r2y = ry - sin(beta_m)*a - sin(beta_m-alpha2)*a;
r1x = rx + cos(beta_m)*a + cos(beta_m+alpha1)*a;
r1y = ry + sin(beta_m)*a + sin(beta_m+alpha1)*a;
switch draw_mode
    case '3d'
        [x, y, z] = ellipsoid(rx,ry,0,a,b,c);
        [x1, y1, z1] = ellipsoid(r2x,r2y,0,a,b,c);
        [x2, y2, z2] = ellipsoid(r1x,r1y,0,a,b,c);
        h_fish.rh = surf(x,y,z,'FaceColor',[0.25 0.59 0.81]);hold on;
        rotate(rh,[0 0 1],rad2deg(beta_m),[rx,ry,0]);
        h_fish.rh1 = surf(x1,y1,z1,'FaceColor',[0.25 0.59 0.81]);
        rotate(rh1,[0 0 1],rad2deg(beta2),[r2x,r2y,0]);
        rh2 = surf(x2,y2,z2,'FaceColor',[0.25 0.59 0.81]);
        h_fish.rotate(rh2,[0 0 1],rad2deg(beta1),[r1x,r1y,0]);
        h_fish.eyeh = plot([r1x+a/2,r1x+a/2],[r1y+b,r1y-b],'ko','MarkerFaceColor',[0.5,0.5,0.5]);
        rotate(h_fish.eyeh,[0 0 1],rad2deg(beta1),[r1x,r1y,0]);
        set(h_fish.rh,'LineStyle','none');
        set(h_fish.rh1,'LineStyle','none');
        set(h_fish.rh2,'LineStyle','none');
        axis equal;
        ax = gca;
        ax.TickLabelInterpreter = 'latex';
        ax.FontSize = 15;
        ax.Box = 'on';
        ax.XColor = 'none';
        ax.YColor = 'none';
        ax.LineWidth = 2;
        view([0 90]);
        camlight left;

    case '2d'

        t = linspace(0,2*pi,100);
        ellipse = [a*cos(t);
            1.25*b*sin(t);];
        rotated = Rot2D(beta_m)*ellipse;
        x = rotated(1,:)+rx;
        y = rotated(2,:)+ry;
        sh = polyshape(x(1:end-1),y(1:end-1));
        rotated = Rot2D(beta1)*ellipse;
        x1 = rotated(1,:)+r1x;
        y1 = rotated(2,:)+r1y;
        sh1 = polyshape(x1(1:end-1),y1(1:end-1));
        rotated = Rot2D(beta2)*ellipse;
        x2 = rotated(1,:)+r2x;
        y2 = rotated(2,:)+r2y;
        sh2 = polyshape(x2(1:end-1),y2(1:end-1));
        eyes = Rot2D(beta1)*[[a/2,a/2];[b,-b]];
        eyes(1,:) = r1x+eyes(1,:);
        eyes(2,:) = r1y+eyes(2,:);
        switch nargin
            case 9
                h_fish.middle = plot(sh,'FaceColor',[0.25 0.59 0.81],'FaceAlpha',0.8,'EdgeColor','none');
                hold on;
                h_fish.head = plot(sh1,'FaceColor',[0.25 0.59 0.81],'FaceAlpha',0.8,'EdgeColor','none');
                h_fish.tail = plot(sh2,'FaceColor',[0.25 0.59 0.81],'FaceAlpha',0.8,'EdgeColor','none');
                h_fish.eyeh = plot(eyes(1,:),eyes(2,:),'ko','MarkerFaceColor',[0.2,0.2,0.2],'MarkerSize',6);
            case 10
                h_fish.middle.Shape = sh;
                h_fish.head.Shape = sh1;
                h_fish.tail.Shape = sh2;
                h_fish.eyeh.XData = eyes(1,:);
                h_fish.eyeh.YData = eyes(2,:);
        end
        drawnow;
end
end
function R = Rot2D(angle)
R = [cos(angle), -sin(angle);
    sin(angle), cos(angle)];
end