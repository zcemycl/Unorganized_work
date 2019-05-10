% Only for Calculation
clear all; close all; clc;

step_y = 0.1; step_x = 0.15; % step_x cannot be greater than 0.2 or == 0.1
d_two_legs = 0.2;
height = 0.0005; L = 0.1;
n = 50; new_ref = 0;
M = 0.1;
RInertia = M*L^2/3;

step_x_step = 0.14;

nfrb = floor(0.5/step_x);

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
% Test force, velocity and displacement
% --------------------------------------------------------------- %
% rear flat
Tlist_1_1 = [];
Tlist_2_1 = [];
Tlist_3_1 = [];
F_1 = [1,0];
for i = 1:n
    Phi = phiIr + phiIMr*i;
    T = tau(Phi, 'r', F_1);
    Tlist_1_1(i) = T(1);
    Tlist_2_1(i) = T(2);
    Tlist_3_1(i) = T(3);
end

Tlist_1_2 = [];
Tlist_2_2 = [];
Tlist_3_2 = [];
F_2 = [-1,0];
for i = 1:n
    Phi = phiMr + phiMFr*i;
    T = tau(Phi, 'r', F_2);
    Tlist_1_2(i) = T(1);
    Tlist_2_2(i) = T(2);
    Tlist_3_2(i) = T(3);
end

k = 3; l = 2;
subplot(l,k,1)
plot(Tlist_1_1); hold on; plot(Tlist_1_2); legend(['1';'2']); hold off
subplot(l,k,2)
plot(Tlist_2_1); hold on; plot(Tlist_2_2); legend(['1';'2']); hold off
subplot(l,k,3)
plot(Tlist_3_1); hold on; plot(Tlist_3_2); legend(['1';'2']); hold off

avlist_1_1 = []; aalist_1_1 = [];
avlist_2_1 = []; aalist_2_1 = [];
avlist_3_1 = []; aalist_3_1 = [];
avlist_1_2 = []; aalist_1_2 = [];
avlist_2_2 = []; aalist_2_2 = [];
avlist_3_2 = []; aalist_3_2 = [];
for i = 1:n
    avlist_1_1(i) = sum(Tlist_1_1(1:i)); 
    avlist_2_1(i) = sum(Tlist_2_1(1:i)); 
    avlist_3_1(i) = sum(Tlist_3_1(1:i)); 
end
for i = 1:n
    avlist_1_2(i) = avlist_1_1(n) + sum(Tlist_1_2(1:i)); 
    avlist_2_2(i) = avlist_2_1(n) + sum(Tlist_2_2(1:i)); 
    avlist_3_2(i) = avlist_3_1(n) + sum(Tlist_3_2(1:i)); 
end
subplot(l,k,4)
plot(avlist_1_1); hold on; plot(avlist_1_2); legend(['1';'2']); hold off
subplot(l,k,5)
plot(avlist_2_1); hold on; plot(avlist_2_2); legend(['1';'2']); hold off
subplot(l,k,6)
plot(avlist_3_1); hold on; plot(avlist_3_2); legend(['1';'2']); hold off
% --------------------------------------------------------------- %                                  
function [PE, MID] = potential(M, Phi, L,mode)
[r1T, r2T, r3T, r4T] = transformation(Phi, L, mode);
mid1 = M(1); mid2 = M(2); mid3 = M(3); mid4 = M(4);
if mode == 'r'
    PE1 = abs((r1T(2))*0.5 - mid1)*0.1*9.81;
    PE2 = abs((r1T(2)+r2T(2))*0.5 - mid2)*0.1*9.81;
    PE3 = abs((r2T(2)+r3T(2))*0.5 - mid3)*0.1*9.81;
    PE4 = abs((r3T(2)+r4T(2))*0.5 - mid4)*0.1*9.81;
elseif mode == 'f'
    PE1 = abs((r1T(2)+r2T(2))*0.5 - mid1)*0.1*9.81;
    PE2 = abs((r2T(2)+r3T(2))*0.5 - mid2)*0.1*9.81;
    PE3 = abs((r4T(2)+r3T(2))*0.5 - mid3)*0.1*9.81;
    PE4 = abs(r4T(2)*0.5 - mid4)*0.1*9.81;
end
MID = midpoint_link(r1T, r2T, r3T, r4T, mode);
PE = PE1+PE2+PE3+PE4;
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

function T = tau(Phi, mode, F)                                   
phi1 = Phi(1); phi2 = Phi(2); phi3 = Phi(3);
L = 0.1; g = 9.81; M = 0.4;
if mode == 'r'
    J_T = L*[cos(phi1)+cos(phi1+phi2)+cos(phi1+phi2+phi3),...
        sin(phi1)+sin(phi1+phi2)+sin(phi1+phi2+phi3);
        cos(phi1+phi2)+cos(phi1+phi2+phi3),...
        sin(phi1+phi2)+sin(phi1+phi2+phi3);
        cos(phi1+phi2+phi3),...
        sin(phi1+phi2+phi3)];
elseif mode == 'f'
    J_T = L*[cos(phi1+phi2+phi3),...
        sin(phi1+phi2+phi3);
        cos(phi2+phi3)+cos(phi1+phi2+phi3),...
        sin(phi2+phi3)+sin(phi1+phi2+phi3);
        cos(phi3)+cos(phi2+phi3)+cos(phi1+phi2+phi3),...
        sin(phi3)+sin(phi2+phi3)+sin(phi1+phi2+phi3)];
end
T = J_T*F';
end

