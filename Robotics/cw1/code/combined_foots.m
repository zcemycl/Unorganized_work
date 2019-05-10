clear all; close all; clc;
% combined foots
step_y = 0.1; n = 25; height = 0.0005;

% initialization
lineObj = animInit();

% --------------------------------------------------- %
% flat surface
% ref_x helps to update the new reference foot location
% ref_y helps to walk on the level above the ground
% --------------------------------------------------- %
% front foot motion
x0f = -0.05; x1f = -0.2;
[Phi1f, Phiintf, Phi0f] = angle(x1f,x0f,step_y,'f');
Phi0intstepf = (Phiintf - Phi0f)/n;
Phiint1stepf = (Phi1f - Phiintf)/n;

% rear foot motion
frontref = 0; % front foot fixed 
x0r = 0.2; x1r = 0.05;
[Phi1r, Phiintr, Phi0r] = angle(x1r,x0r,step_y,'r');
Phi0intstepr = (Phiintr - Phi0r)/n;
Phiint1stepr = (Phi1r - Phiintr)/n;
% --------------------------------------------------- %
% staircase motion
% --------------------------------------------------- %
% front foot motion
x1 = -0.15; y1 = 0.0005;
x0 = -0.05; xintstep = 0.5*(x0+x1); 
yintstep = y1 + step_y;

[phi1stepf1, phi2stepf1, phi3stepf1] = pure_angle(x1,y1,'f');
Phistepf1 = [phi1stepf1, phi2stepf1, phi3stepf1];

[phi1stepf0, phi2stepf0, phi3stepf0] = pure_angle(x0,0,'f');
Phistepf0 = [phi1stepf0, phi2stepf0, phi3stepf0];

[phi1stepfint, phi2stepfint, phi3stepfint] = pure_angle(xintstep,yintstep,'f');
Phistepfint = [phi1stepfint, phi2stepfint, phi3stepfint];

Phisf0fint = (Phistepfint - Phistepf0)/n;
Phifintsf1 = (Phistepf1 - Phistepfint)/n;

% rear foot motion
x0stepr = 0.15; y0stepr = -0.0005;
x1stepr = 0.03; y1stepr = 0;
xintstepr = 0.5*(x0stepr + x1stepr);
yintstepr = y0stepr + step_y;
[phi1stepr1, phi2stepr1, phi3stepr1] = pure_angle(x1stepr,y1stepr,'r');
Phistepr1 = [phi1stepr1, phi2stepr1, phi3stepr1];

[phi1stepr0, phi2stepr0, phi3stepr0] = pure_angle(x0stepr,y0stepr,'r');
Phistepr0 = [phi1stepr0, phi2stepr0, phi3stepr0];

[phi1steprint, phi2steprint, phi3steprint] = pure_angle(xintstepr,yintstepr,'r');
Phisteprint = [phi1steprint, phi2steprint, phi3steprint];

Phisf0rint = (Phisteprint - Phistepr0)/n;
Phirintsf1 = (Phistepr1 - Phisteprint)/n;
% --------------------------------------------------- %
% Animation part
animation(Phi0r, Phi0intstepr, frontref, 0, n, lineObj, 'r');
new_ref = animation(Phiintr, Phiint1stepr, frontref, 0, n, lineObj, 'r');
for j = 1:3
animation(Phi0f, Phi0intstepf, new_ref, 0, n, lineObj, 'f');
new_ref = animation(Phiintf, Phiint1stepf, new_ref, 0, n, lineObj, 'f');
animation(Phi0r, Phi0intstepr, new_ref, 0, n, lineObj, 'r');
new_ref = animation(Phiintr, Phiint1stepr, new_ref, 0, n, lineObj, 'r');
end
animation(Phistepf0, Phisf0fint, new_ref, 0, n, lineObj, 'f');
new_ref = animation(Phistepfint, Phifintsf1, new_ref, 0, n, lineObj, 'f');
animation(Phistepr0, Phisf0rint, new_ref, height, n, lineObj, 'r');
new_ref = animation(Phisteprint, Phirintsf1, new_ref, height, n, lineObj, 'r');
for j = 1:4
animation(Phi0f, Phi0intstepf, new_ref, height, n, lineObj, 'f');
new_ref = animation(Phiintf, Phiint1stepf, new_ref, height, n, lineObj, 'f');
animation(Phi0r, Phi0intstepr, new_ref, height, n, lineObj, 'r');
new_ref = animation(Phiintr, Phiint1stepr, new_ref, height, n, lineObj, 'r');
end

% -------------------------------------------------------- %
% Energy calculation
% Energy due to gravity
% rear foot flat surface energy
yCM_prev = yCMcal(Phi0r);
Egtotal1 = 0;
for i = 1:n
    Phi = Phi0r +  Phi0intstepr*i;      
    [yCM, Eg] = gravityenergy(yCM_prev, Phi);
    yCM_prev = yCM;
    Egtotal1 = Egtotal1 + Eg;
end
Egtotal2 = 0; 
yCM_prev2 = yCMcal(Phiintr);
for i = 1:n
    Phi = Phiintr + Phiint1stepr*i;    
    [yCM, Eg] = gravityenergy(yCM_prev2, Phi);
    yCM_prev2 = yCM;
    Egtotal2 = Egtotal2 + Eg;
end
% front foot flat surface energy
Egtotal3 = 0;
yCM_prev3 = yCMcal(Phi0f);
for i = 1:n
    Phi = Phi0f + Phi0intstepf*i;      
    [yCM, Eg] = gravityenergy(yCM_prev3, Phi);
    yCM_prev3 = yCM;
    Egtotal3 = Egtotal3 + Eg;
end
Egtotal4 = 0;
yCM_prev4 = yCMcal(Phiintf);
for i = 1:n
    Phi = Phiintf + Phiint1stepf*i;      
    [yCM, Eg] = gravityenergy(yCM_prev4, Phi);
    yCM_prev4 = yCM;
    Egtotal4 = Egtotal4 + Eg;
end
% Total energy due to motors
% -------------------------------------------------------- %
totalenergy1 = 0;
for i = 1:n
Phi = Phi0r + Phi0intstepr*i;
total = totalcal(Phi, 'r', Phi0intstepr);
totalenergy1 = totalenergy1 + total;
end
% -------------------------------------------------------- %
% Functions
function lineObj = animInit()
hold on
axis equal
axis([-1.2 0.3 0 0.25])
plot([-1.2 -0.5 -0.5],[0.0005 0.0005 0],...
         'k','LineWidth',3);


%link line objects
lineObj.h1T = line(0,0,'color','k','LineWidth',2);
lineObj.h2T = line(0,0,'color','k','LineWidth',2);
lineObj.h3T = line(0,0,'color','k','LineWidth',2);
lineObj.h4T = line(0,0,'color','k','LineWidth',2);
end 

function [new_ref] = animation(Phi0, Phi0intstep, x_ref, y_ref, n, lineObj, mode)
L = 0.1;

for i = 1:n
    if mode == 'r'
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
        set(lineObj.h1T,'xdata',[x_ref x_ref+r1T(1)],'ydata',[y_ref y_ref+r1T(2)])
        set(lineObj.h2T,'xdata',[x_ref+r1T(1) x_ref+r2T(1)],'ydata',[y_ref+r1T(2) y_ref+r2T(2)])
        set(lineObj.h3T,'xdata',[x_ref+r2T(1) x_ref+r3T(1)],'ydata',[y_ref+r2T(2) y_ref+r3T(2)])
        set(lineObj.h4T,'xdata',[x_ref+r3T(1) x_ref+r4T(1)],'ydata',[y_ref+r3T(2) y_ref+r4T(2)])
        drawnow;
        
        new_ref = x_ref+r4T(1);
        
    elseif mode == 'f'
        phi1_T = Phi0(1)+Phi0intstep(1)*i;
        phi2_T = Phi0(2)+Phi0intstep(2)*i;
        phi3_T = Phi0(3)+Phi0intstep(3)*i;
        A10T = [1 , 0 ,  0; ...
               0 , 1 , L; ...
               0 , 0 , 1 ];

        A21T = [cos(phi1_T) , -sin(phi1_T) , -L*sin(phi1_T) ; ...
                sin(phi1_T) ,  cos(phi1_T) , L*cos(phi1_T); ...
                0        , 0          , 1  ];

        A32T = [cos(phi2_T) , -sin(phi2_T) , -L*sin(phi2_T); ...
                sin(phi2_T) ,  cos(phi2_T) , L*cos(phi2_T); ...
               0         , 0         , 1 ];
        A43T = [cos(phi3_T) , -sin(phi3_T) , -L*sin(phi3_T); ...
                sin(phi3_T) ,  cos(phi3_T) , L*cos(phi3_T); ...
                0         , 0         , 1 ];
        r1T = A10T*A43T*A32T*A21T;   %position link 1
        r1T = r1T(:,3);
        r1T = r1T(1:2);
        r2T = A10T*A43T*A32T;   %position link 2
        r2T = r2T(:,3);
        r2T = r2T(1:2);
        r3T = A10T*A43T;
        r3T = r3T(:,3);
        r3T = r3T(1:2);
        r4T = A10T;
        r4T = r4T(:,3);
        r4T = r4T(1:2);

        %set target position and pause
        set(lineObj.h1T,'xdata',[x_ref+r1T(1) x_ref+r2T(1)],'ydata',[y_ref+r1T(2) y_ref+r2T(2)])
        set(lineObj.h2T,'xdata',[x_ref+r2T(1) x_ref+r3T(1)],'ydata',[y_ref+r2T(2) y_ref+r3T(2)])
        set(lineObj.h3T,'xdata',[x_ref+r4T(1) x_ref+r3T(1)],'ydata',[y_ref+r4T(2) y_ref+r3T(2)])
        set(lineObj.h4T,'xdata',[x_ref x_ref+r4T(1)],'ydata',[y_ref y_ref+r4T(2)])
        drawnow;
        
        new_ref = x_ref+r1T(1);
    
    end
    
end

end

function [Phi1, Phiint, Phi0]= angle(x1, x0, step_y, mode)
xint = (x0 + x1)/2;
[Phi01, Phi02, Phi03] = pure_angle(x0,0, mode);
[Phi11, Phi12, Phi13] = pure_angle(x1,0, mode);
[Phiint1, Phiint2, Phiint3] = pure_angle(xint,step_y, mode);

Phi0 = [Phi01, Phi02, Phi03];
Phi1 = [Phi11, Phi12, Phi13];
Phiint = [Phiint1, Phiint2, Phiint3];
end

function [phi1, phi2, phi3] = pure_angle(x,y,mode)
L = 0.1; 
phi2 = acos((x^2 + y^2)/2/L^2 - 1);
if mode == 'f'
    phi1 = asin(0.5*(-x/L + sin(phi2)*y/L/(1+cos(phi2))));
    phi3 = pi - phi1 - phi2;
elseif mode == 'r'
    phi3 = asin(0.5*(x/L + sin(phi2)*((y-L)/L +1)/(cos(phi2)+1)));
    phi1 = pi - phi2 - phi3; 
end

end 



function [yCM, Eg] = gravityenergy(yCM_prev, Phi)
    L = 0.1; g = 9.81; M = 0.4;
    phi1 = Phi(1); phi2 = Phi(2); phi3 = Phi(3);
    % centre of mass along y-axis
    yCM = L/8*(7+5*cos(phi1)+3*cos(phi1+phi2)+cos(phi1+phi2+phi3));
    Eg = abs(yCM - yCM_prev)*M*g;
end

function yCM = yCMcal(Phi)
L = 0.1; g = 9.81; M = 0.4;
phi1 = Phi(1); phi2 = Phi(2); phi3 = Phi(3);
yCM = L/8*(7+5*cos(phi1)+3*cos(phi1+phi2)+cos(phi1+phi2+phi3));
end
function total = totalcal(Phi, mode, dtheta)
phi1 = Phi(1); phi2 = Phi(2); phi3 = Phi(3);
L = 0.1; g = 9.81; M = 0.4;
if mode == 'r'
    tau1 = M*g*(sin(phi1)+sin(phi1+phi2)+sin(phi1+phi2+phi3));
    tau2 = M*g*(sin(phi1+phi2)+sin(phi1+phi2+phi3));
    tau3 = M*g*(sin(phi1+phi2+phi3));
elseif mode == 'f'
    tau1 = M*g*(sin(phi1+phi2+phi3));
    tau2 = M*g*(sin(phi2+phi3)+sin(phi1+phi2+phi3));
    tau3 = M*g*(sin(phi3)+sin(phi2+phi3)+sin(phi1+phi2+phi3));
end
torquevec = abs([tau1, tau2, tau3]);
total = torquevec*abs(dtheta)';
end 
% function dispangle(Phi)
% disp(Phi*180/pi);
% end

