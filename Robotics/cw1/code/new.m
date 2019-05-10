% new
close all; clc; clear all;
x0 = 0.2; x1 = 0.05;
step_y = 0.025;n = 150;
[Phi1, Phiint, Phi0] = angle(x1,x0,step_y);
dispangle(Phi1);
Phi0intstep = (Phiint - Phi0)/n;
Phiint1step = (Phi1 - Phiint)/n;
lineObj = animInit();
animation(Phi0, Phi0intstep, n, lineObj)
animation(Phiint, Phiint1step, n, lineObj)
function lineObj = animInit()
hold on
axis equal
axis([-0.4 0.4 0 0.45])

%link line objects
lineObj.h1T = line(0,0,'color','k','LineWidth',2);
lineObj.h2T = line(0,0,'color','k','LineWidth',2);
lineObj.h3T = line(0,0,'color','k','LineWidth',2);
lineObj.h4T = line(0,0,'color','k','LineWidth',2);
%joint line objects
lineObj.d1T = line(0,0,'color','k','LineWidth',5,'Marker','o');
lineObj.d2T = line(0,0,'color','k','LineWidth',5,'Marker','o');
lineObj.d3T = line(0,0,'color','k','LineWidth',5,'Marker','o');
lineObj.d4T = line(0,0,'color','k','LineWidth',5,'Marker','o');
end 
function animation(Phi0, Phi0intstep, n, lineObj)
L = 0.1;
for i = 1:n
    phi1_T = Phi0(1)+Phi0intstep(1)*i;
    phi2_T = Phi0(2)+Phi0intstep(2)*i;
    phi3_T = Phi0(3)+Phi0intstep(3)*i;
    A10T = [1 , 0 ,  0; ...
           0 , 1 , L; ...
           0 , 0 , 1 ];

    A21T = [cos(phi1_T) , sin(phi1_T) , L*sin(phi1_T) ; ...
           -sin(phi1_T) ,  cos(phi1_T) , L*cos(phi1_T); ...
            0        , 0          , 1  ];

    A32T = [cos(phi2_T) , sin(phi2_T) , L*sin(phi2_T); ...
           -sin(phi2_T) ,  cos(phi2_T) , L*cos(phi2_T); ...
           0         , 0         , 1 ];
    A43T = [cos(phi3_T) , sin(phi3_T) , L*sin(phi3_T); ...
           -sin(phi3_T) ,  cos(phi3_T) , L*cos(phi3_T); ...
            0         , 0         , 1 ];
    r1T = A10T(:,3);   %position link 1
    r1T = r1T(1:2);
    r2T = A10T*A21T;   %position link 2
    r2T = r2T(:,3);
    r2T = r2T(1:2);
    r3T = A10T*A21T*A32T;
    r3T = r3T(:,3);
    r3T = r3T(1:2);
    r4T = A10T*A21T*A32T*A43T;
    r4T = r4T(:,3);
    r4T = r4T(1:2);
    
    
    %set target position and pause
    set(lineObj.h1T,'xdata',[0 r1T(1)],'ydata',[0 r1T(2)])
    set(lineObj.h2T,'xdata',[r1T(1) r2T(1)],'ydata',[r1T(2) r2T(2)])
    set(lineObj.h3T,'xdata',[r2T(1) r3T(1)],'ydata',[r2T(2) r3T(2)])
    set(lineObj.h4T,'xdata',[r3T(1) r4T(1)],'ydata',[r3T(2) r4T(2)])
    set(lineObj.d1T,'xdata',r1T(1),'ydata',r1T(2))
    set(lineObj.d2T,'xdata',r2T(1),'ydata',r2T(2))
    set(lineObj.d3T,'xdata',r3T(1),'ydata',r3T(2))
    drawnow;
end
end

function [Phi1, Phiint, Phi0]= angle(x1, x0, step_y)
xint = (x0 + x1)/2;
[Phi01, Phi02, Phi03] = pure_angle(x0,0);
Phi0 = [Phi01, Phi02, Phi03];
[Phi11, Phi12, Phi13] = pure_angle(x1,0);
Phi1 = [Phi11, Phi12, Phi13];
[Phiint1, Phiint2, Phiint3] = pure_angle(xint,step_y);
Phiint = [Phiint1, Phiint2, Phiint3];
end

function [phi1, phi2, phi3] = pure_angle(x,y)
L = 0.1;
phi2 = acos((x^2+(y-L)^2)/(2*L^2) + (y-L)/L -0.5);
phi3 = asin(0.5*(x/L + sin(phi2)*((y-L)/L +1)/(cos(phi2)+1)));
phi1 = pi - phi2 - phi3;
end

function dispangle(Phi)
disp(Phi*180/pi);
end