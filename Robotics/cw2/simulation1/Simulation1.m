clear all; close all; clc;

step_y = 0.08; step_x = 0.15; % step_x cannot be greater than 0.2 or == 0.1
d_two_legs = 0.2;
height = 0.05; L = 0.14;
Lf = 0.2;
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
