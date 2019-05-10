clear all;close all; clc;
% Settings
% Predict joint angles trajectory --> Set time required
% --> torque required --> Energy efficiency
% Constraint
% - Sum of the joint angles = 180
% Trajectory Profie
% angle = angle_acceleration*t_c*(t_f-t_c)

Tsample = 0.1; Taction = 20;
n = Taction/Tsample;

x_sep = 0.2; x_step = 0.15;
y_step = 0.05; height = 0.0005;
L = 0.1; M = 0.1; I  = M*L^2/3;
new_ref = 0;

x_stair = sqrt(x_sep^2 - (0.5*(height+2*y_step))^2);
x_stair_floor = floor(x_stair);

a_index = 0.5;

nf = floor(0.5/x_step); nr = floor((x_sep+0.5)/x_step);

% initialization
lineObj = lim(x_sep, x_step, y_step);
% --------------------------------------------------------------- %
% Motion calculation (Trajectory)
% --------------------------------------------------------------- %
% Rear Flat
xirf = x_sep; xfrf = x_sep-x_step;
yirf = 0; yfrf = 0;
[phiIrf, phiMrf, phiFrf, phiIMrf, phiMFrf] = angle_step(xirf, ...
    xfrf, yirf, yfrf, y_step, n, 'r');
% Front Flat
xiff = -xfrf; xfff = -x_sep;
yiff = 0; yfff = 0;
[phiIff, phiMff, phiFff, phiIMff, phiMFff] = angle_step(xiff, ...
    xfff, yiff, yfff, y_step, n, 'f');
% Front Step
xifs = -xfrf; xffs = -0.195;
yifs = 0; yffs = height;
[phiIfs, phiMfs, phiFfs, phiIMfs, phiMFfs] = angle_step(xifs, ...
    xffs, yifs, yffs, y_step, n, 'f');

% --------------------------------------------------------------- %
% Animation: Flat 
% --------------------------------------------------------------- %
for i = 1:nf
[T1r,T2r,T3r,new_ref,new_ref0] = motion(phiIrf, phiMrf, phiIMrf, phiMFrf, 'r',...
    L, lineObj, new_ref, 0, n, a_index, Taction, I);
[T1f,T2f,T3f,new_ref,new_ref0] = motion(phiIff, phiMff, phiIMff, phiMFff, 'f',...
    L, lineObj, new_ref, 0, n, a_index, Taction, I);
end
% --------------------------------------------------------------- %
% Animation: Flat 
% --------------------------------------------------------------- %
if new_ref < new_ref0-x_step
    [T1r,T2r,T3r,new_ref,new_ref0] = motion(phiIrf, phiMrf, phiIMrf, phiMFrf, 'r',...
            L, lineObj, new_ref, 0, n, a_index, Taction, I);
end
if new_ref-(-0.5) >= 0.14
    [T1f,T2f,T3f,new_ref,new_ref0] = motion(phiIff, phiMff, phiIMff, phiMFff, 'f',...
    L, lineObj, new_ref, 0, n, a_index, Taction, I);
end
    
% Rare Step
xirs = -0.36+0.555;  xfrs = (-0.5 - (-0.555))/2;
yirs = -height; yfrs = 0;
[phiIrs, phiMrs, phiFrs, phiIMrs, phiMFrs] = angle_step(xirs, ...
    xfrs, yirs, yfrs, y_step, n, 'r');
% Front Adjust
xifa = -0.0275; xffa = x_step-x_sep;
yifa = 0; yffa = 0;
[phiIfa, phiMfa, phiFfa, phiIMfa, phiMFfa] = angle_step(xifa, ...
    xffa, yifa, yffa, y_step, n, 'f');
[T1fs,T2fs,T3fs,new_ref] = motion(phiIfs, phiMfs, phiIMfs, phiMFfs, 'f',...
    L, lineObj, new_ref, 0, n, a_index, Taction, I);
[T1rs,T2rs,T3rs,new_ref] = motion(phiIrs, phiMrs, phiIMrs, phiMFrs, 'r',...
    L, lineObj, new_ref, height, n, a_index, Taction, I);
new_ref_2 = new_ref;
[T1fa,T2fa,T3fa,new_ref] = motion(phiIfa, phiMfa, phiIMfa, phiMFfa, 'f',...
    L, lineObj, new_ref, height, n, a_index, Taction, I);
[T1f,T2f,T3f,new_ref] = motion(phiIff, phiMff, phiIMff, phiMFff, 'f',...
    L, lineObj, new_ref_2, height, n, a_index, Taction, I);
for i = 1:nf
[T1r,T2r,T3r,new_ref] = motion(phiIrf, phiMrf, phiIMrf, phiMFrf, 'r',...
    L, lineObj, new_ref, 0, n, a_index, Taction, I);
[T1f,T2f,T3f,new_ref] = motion(phiIff, phiMff, phiIMff, phiMFff, 'f',...
    L, lineObj, new_ref, 0, n, a_index, Taction, I);
end


% --------------------------------------------------------------- %
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

function T = torquecal(a_index, t, I, Phi)
    T = I/t^2/a_index/(1-a_index)*Phi;
end

function lineObj = lim(x_sep, x_step, y_step)
x_int = (x_sep+x_step)/2;
r = x_int^2 + y_step^2;
if r <= x_sep^2
    lineObj = animInit();
end
end

function [T1,T2,T3,new_ref, xref0] = motion(phiI, phiM, phiIM, ...
    phiMF, mode, L,...
    lineObj, xref, yref, n, a_index, t, I)
T1 = []; T2 = []; T3 = [];
if mode == 'r'
    T = torquecal(a_index, t, I, phiM-phiI);
    T1(1) = T(1); T2(1) = T(2); T3(1) = T(3);
    for i = 1:n
        Phi = phiI + phiIM*i;
        [r1T, r2T, r3T, r4T] = transformation(Phi, L, 'r');
        animation(r1T, r2T, r3T, r4T, xref, yref, lineObj, 'r');
    end
    T = torquecal(a_index, t, I, phiMF*n);
    T1(2) = T(1); T2(2) = T(2); T3(2) = T(3);
    for i = 1:n
        Phi2 = phiM + phiMF*i;
        [r1T, r2T, r3T, r4T] = transformation(Phi2, L, 'r');
        [new_ref] = animation(r1T, r2T, r3T, r4T, xref, yref, lineObj, 'r');
    end
elseif mode == 'f'
    T = torquecal(a_index, t, I, phiM-phiI);
    T1(1) = T(1); T2(1) = T(2); T3(1) = T(3);
    for i = 1:n
        Phi = phiI + phiIM*i;
        [r1T, r2T, r3T, r4T] = transformation(Phi, L, 'f');
        animation(r1T, r2T, r3T, r4T, xref, yref, lineObj, 'f');
    end
    T = torquecal(a_index, t, I, phiMF*n);
    T1(2) = T(1); T2(2) = T(2); T3(2) = T(3);
    for i = 1:n
        Phi2 = phiM + phiMF*i;
        [r1T, r2T, r3T, r4T] = transformation(Phi2, L, 'f');
        [new_ref] = animation(r1T, r2T, r3T, r4T, xref, yref, lineObj, 'f');
    end
end
xref0 = xref;
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


