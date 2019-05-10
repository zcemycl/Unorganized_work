clear all; close all; clc;

step_y = 0.08; step_x = 0.15; % step_x cannot be greater than 0.2 or == 0.1
d_two_legs = 0.2;
height = 0.05; L = 0.14;
L0 = 0.2;
L1 = 0.14;
L2 = 0.13;
L3 = 0.15;
L4 = 0.13;
L5 = 0.2;
n = 20; new_ref = 0;
M = 0.1;
RInertia = M*L^2/3;

step_x_step = 0.14;

nfrb = floor(0.5/step_x);


% initialization
lineObj = animInit(height);

% --------------------------------------------------------------- %
% Motion trajectory
% --------------------------------------------------------------- %
% rear flat
xir = d_two_legs; xfr = xir - step_x; 

[phiIr, phiMr, phiFr, phiIMr, phiMFr] = angle_step(xir, xfr,...
                                        0, 0, step_y, n, 'r');
% front flat
xif = -xfr; xff = xif - step_x; 
[phiIf, phiMf, phiFf, phiIMf, phiMFf] = angle_step(xif, xff,...
                                        0, 0, step_y, n, 'f');

% front step                                    
xifs = step_x - d_two_legs; xffs = xifs - step_x_step; 
[phiIfs, phiMfs, phiFfs, phiIMfs, phiMFfs] = angle_step(xifs, xffs,...
                                        0, height, step_y, n, 'f');    
% rear step
xirs = -xffs; xfrs = 0.05; 
[phiIrs, phiMrs, phiFrs, phiIMrs, phiMFrs] = angle_step(xirs, xfrs,...
                                        -height, 0, step_y, n, 'r');

% --------------------------------------------------------------- %
% Motion Animation
% --------------------------------------------------------------- %

                                    
for j = 1:nfrb
new_ref = motion(phiIr, phiMr, phiIMr, phiMFr, 'r', L, lineObj, new_ref, 0, n);
new_ref = motion(phiIf, phiMf, phiIMf, phiMFf, 'f', L, lineObj, new_ref, 0, n);
end
new_ref = motion(phiIr, phiMr, phiIMr, phiMFr, 'r', L, lineObj, new_ref, 0, n);
new_ref = motion(phiIfs, phiMfs, phiIMfs, phiMFfs, 'f', L, lineObj, new_ref, 0, n);
new_ref = motion(phiIrs, phiMrs, phiIMrs, phiMFrs, 'r', L, lineObj, new_ref, height, n);
for k = 1:3
new_ref = motion(phiIf, phiMf, phiIMf, phiMFf, 'f', L, lineObj, new_ref, height, n);
new_ref = motion(phiIr, phiMr, phiIMr, phiMFr, 'r', L, lineObj, new_ref, height, n);
end
new_ref = motion(phiIf, phiMf, phiIMf, phiMFf, 'f', L, lineObj, new_ref, height, n);



% ------------------------------------------------5-------- %
% Functions
% --------------------------------------------------------------- %
function lineObj = animInit(height)
hold on
axis equal
axis([-1.2 0.3 0 0.25])
plot([-1.2 -0.5 -0.5],[height height 0],...
         'k','LineWidth',3);

%link line objects
lineObj.h1T = line(0,0,'color','k','LineWidth',2);
lineObj.h2T = line(0,0,'color','k','LineWidth',2);
lineObj.h3T = line(0,0,'color','k','LineWidth',2);
lineObj.h4T = line(0,0,'color','k','LineWidth',2);
end 

function new_ref = motion(phiI, phiM, phiIM, phiMF, mode, L, lineObj, xref, yref, n)
if mode == 'r'
    for i = 1:n
        Phi = phiI + phiIM*i;
        [r1T, r2T, r3T, r4T] = transformation(Phi, L, 'r');
        animation(r1T, r2T, r3T, r4T, xref, yref, lineObj, 'r');
    end
    
    for i = 1:n
        Phi2 = phiM + phiMF*i;
        [r1T, r2T, r3T, r4T] = transformation(Phi2, L, 'r');
        [new_ref] = animation(r1T, r2T, r3T, r4T, xref, yref, lineObj, 'r');
    end
elseif mode == 'f'
    for i = 1:n
        Phi = phiI + phiIM*i;
        [r1T, r2T, r3T, r4T] = transformation(Phi, L, 'f');
        animation(r1T, r2T, r3T, r4T, xref, yref, lineObj, 'f');
    end
    
    for i = 1:n
        Phi2 = phiM + phiMF*i;
        [r1T, r2T, r3T, r4T] = transformation(Phi2, L, 'f');
        [new_ref] = animation(r1T, r2T, r3T, r4T, xref, yref, lineObj, 'f');
    end
end
end

function [r1T, r2T, r3T, r4T] = transformation(Phi, L, mode)
phi1_T = Phi(1); phi2_T = Phi(2); phi3_T = Phi(3);
A10T = [1 , 0 ,  0; ...
        0 , 1 , L; ...
        0 , 0 , 1 ];
if mode == 'r'
    A21T = [cos(phi1_T) , sin(phi1_T) , L*sin(phi1_T) ; ...
           -sin(phi1_T) ,  cos(phi1_T) , L*cos(phi1_T); ...
            0        , 0          , 1  ];

    A32T = [cos(phi2_T) , sin(phi2_T) , L*sin(phi2_T); ...
           -sin(phi2_T) ,  cos(phi2_T) , L*cos(phi2_T); ...
           0         , 0         , 1 ];
    A43T = [cos(phi3_T) , sin(phi3_T) , L*sin(phi3_T); ...
           -sin(phi3_T) ,  cos(phi3_T) , L*cos(phi3_T); ...
            0         , 0         , 1 ];
elseif mode == 'f'
    A21T = [cos(phi1_T) , -sin(phi1_T) , -L*sin(phi1_T) ; ...
            sin(phi1_T) ,  cos(phi1_T) , L*cos(phi1_T); ...
            0        , 0          , 1  ];

    A32T = [cos(phi2_T) , -sin(phi2_T) , -L*sin(phi2_T); ...
            sin(phi2_T) ,  cos(phi2_T) , L*cos(phi2_T); ...
            0         , 0         , 1 ];
    A43T = [cos(phi3_T) , -sin(phi3_T) , -L*sin(phi3_T); ...
            sin(phi3_T) ,  cos(phi3_T) , L*cos(phi3_T); ...
            0         , 0         , 1 ];
end
if mode == 'r'
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
elseif mode == 'f'
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
end
end

function [new_ref] = animation(r1T, r2T, r3T, r4T, x_ref, y_ref, lineObj, mode)
if mode == 'r'
    %set target position and pause
    set(lineObj.h1T,'xdata',[x_ref x_ref+r1T(1)],'ydata',[y_ref y_ref+r1T(2)])
    set(lineObj.h2T,'xdata',[x_ref+r1T(1) x_ref+r2T(1)],'ydata',[y_ref+r1T(2) y_ref+r2T(2)])
    set(lineObj.h3T,'xdata',[x_ref+r2T(1) x_ref+r3T(1)],'ydata',[y_ref+r2T(2) y_ref+r3T(2)])
    set(lineObj.h4T,'xdata',[x_ref+r3T(1) x_ref+r4T(1)],'ydata',[y_ref+r3T(2) y_ref+r4T(2)])
    drawnow;

    new_ref = x_ref+r4T(1);
elseif mode == 'f'
    %set target position and pause
    set(lineObj.h1T,'xdata',[x_ref+r1T(1) x_ref+r2T(1)],'ydata',[y_ref+r1T(2) y_ref+r2T(2)])
    set(lineObj.h2T,'xdata',[x_ref+r2T(1) x_ref+r3T(1)],'ydata',[y_ref+r2T(2) y_ref+r3T(2)])
    set(lineObj.h3T,'xdata',[x_ref+r4T(1) x_ref+r3T(1)],'ydata',[y_ref+r4T(2) y_ref+r3T(2)])
    set(lineObj.h4T,'xdata',[x_ref x_ref+r4T(1)],'ydata',[y_ref y_ref+r4T(2)])
    drawnow;

    new_ref = x_ref+r1T(1);
end
end

function [phiI, phiM, phiF, phiIM, phiMF] = angle_step(xi, xf, yi, yf, step_y, n, mode)
    phiI = pure_angle(xi, yi, mode);
    phiM = pure_angle(0.5*(xi+xf), 0.5*(yi+yf+2*step_y), mode);
    phiF = pure_angle(xf, yf, mode);
    
    phiIM = (phiM - phiI)/n;
    phiMF = (phiF - phiM)/n;
end

function Phi = pure_angle(x,y,mode)
L = 0.1; 
phi2 = acos((x^2 + y^2)/2/L^2 - 1);
if mode == 'f'
    phi1 = asin(0.5*(-x/L + sin(phi2)*y/L/(1+cos(phi2))));
    phi3 = pi - phi1 - phi2;
elseif mode == 'r'
    phi3 = asin(0.5*(x/L + sin(phi2)*((y-L)/L +1)/(cos(phi2)+1)));
    phi1 = pi - phi2 - phi3; 
end
Phi = [phi1, phi2, phi3];

end 

function MID = midpoint_link(r1T, r2T, r3T, r4T, mode)
if mode == 'r'
    mid1 = (r1T(2))*0.5;
    mid2 = (r1T(2)+r2T(2))*0.5;
    mid3 = (r2T(2)+r3T(2))*0.5;
    mid4 = (r3T(2)+r4T(2))*0.5;
elseif mode == 'f'
    mid1 = (r1T(2)+r2T(2))*0.5;
    mid2 = (r2T(2)+r3T(2))*0.5;
    mid3 = (r4T(2)+r3T(2))*0.5;
    mid4 = r4T(2)*0.5;
end
MID = [mid1, mid2, mid3, mid4];
end
