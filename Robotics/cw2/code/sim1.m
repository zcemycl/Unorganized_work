%% Simulation
clear all; close all; clc;
%% Constant
% Time frame
T = 0.1;
% gravity
g = 10;
% Robot info.
% structure
L_foot  = 0.2 ; L_leg   = 0.14;
front_x = 0   ; front_y = 0   ;
rear_x  = 0.28; rear_y  = 0   ;
max_sep = 0.28; min_sep = 0.2 ;
% motion
max_x   = -0.08; 
% min dist between robot and platform
d       = 0.02;
% Torque
Tor = 15;

% Motor 
s = 0.04;
m = 0.1;
% Platform info
plat_sx = -0.5; 
plat_ex = -1.5;
plat_h  = 0.05;

%% Auto - Calculation
% Total distance before the platform
% front
df_bp = plat_sx + d - front_x;
num   = df_bp/max_x          ; % num of steps from start to plat

% All position walk 1 on flat 
movecoor1 = position(front_x,rear_x, plat_h, ...
    max_x, num);
polar = pos_angle(movecoor1, L_leg);
expandpolar = expandang(polar);
Rpos  = transformation(expandpolar, L_leg);

%% Animation
lineObj = animInit(plat_sx, plat_ex, plat_h);
animation(Rpos, movecoor1, expandpolar, lineObj);
%% Functions
function lineObj = animInit(plat_sx, plat_ex, plat_h)
hold on
axis equal
axis([-2 1 0 0.35])
plot([plat_ex plat_ex plat_sx plat_sx],[0 plat_h plat_h 0],...
         'k','LineWidth',3);

%link line objects
lineObj.h0T = line(0,0,'color','k','LineWidth',2);
lineObj.h1T = line(0,0,'color','k','LineWidth',2);
lineObj.h2T = line(0,0,'color','k','LineWidth',2);
lineObj.h3T = line(0,0,'color','k','LineWidth',2);
lineObj.h4T = line(0,0,'color','k','LineWidth',2);
lineObj.h5T = line(0,0,'color','k','LineWidth',2);
end 
function movecoor = position(frontx,rearx, height, distance, num)
movecoor = zeros(4,4*num+3);
movecoor(1,1) = frontx;
movecoor(2,1) = rearx;
dim = size(movecoor);
for i = 2:dim(2)
    if rem(i,2) == 0
        if rem(i,4) == 0
            % Front %
            % arrange x
            movecoor(1,i) = movecoor(1,i-1)+0.5*distance;
            movecoor(1,i+1:i+3) = movecoor(1,i-1)+distance;
            % arrange y
            movecoor(3,i) = height;
        else
            % Rear %
            % arrange x
            movecoor(2,i) = movecoor(2,i-1)+0.5*distance;
            movecoor(2,i+1:i+3) = movecoor(2,i-1)+distance;
            % arrange y
            movecoor(4,i) = height;
        end
    end
end
movecoor = movecoor(:,1:dim(2));
end

function polar = pos_angle(coor, L)
x = coor(2,:)-coor(1,:);
y = coor(4,:)-coor(3,:);
phi2 = acos((x.^2 + y.^2)./2./L^2 - 1);
phi3 = asin(0.5.*(x./L + sin(phi2).*((y-L)./L +1)./(cos(phi2)+1)));
phi1 = pi - phi2 - phi3; 
polar = [phi1;phi2;phi3];

dim = size(polar);
leg = zeros(1,dim(2));
for i = 1:dim(2)
    if rem(i,2) == 0
        if rem(i,4) == 0
            for k = 0:1
                leg(i+k) = 'f';% 102 (ASCII)
            end
        else
            for k = 0:1
                leg(i+k) = 'r';% 114 (ASCII)         
            end
        end
    end
end
leg(1) = 'r';
polar = [polar;leg];

end
%----------------------------------------------------- %
function R_all = transformation(ang, L)
dim = size(ang);
R_all = [];
for i = 1:dim(2)
    phi1_T = ang(1,i); phi2_T = ang(2,i); phi3_T = ang(3,i);
    state = ang(4,i);
    if state == 114
        A10T = [1 , 0 ,  0; ...
                0 , 1 , L; ...
                0 , 0 , 1 ];
        A21T = [cos(phi1_T) , sin(phi1_T) , L*sin(phi1_T) ; ...
               -sin(phi1_T) ,  cos(phi1_T) , L*cos(phi1_T); ...
                0        , 0          , 1  ];

        A32T = [cos(phi2_T) , sin(phi2_T) , L*sin(phi2_T); ...
               -sin(phi2_T) ,  cos(phi2_T) , L*cos(phi2_T); ...
               0          , 0         , 1 ];
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
    elseif state == 102
        A34T = [1 , 0 ,  0; ...
                0 , 1 , L; ...
                0 , 0 , 1 ];
        A23T = [cos(phi3_T) , -sin(phi3_T) , -L*sin(phi3_T); ...
                sin(phi3_T) ,  cos(phi3_T) , L*cos(phi3_T); ...
                0         , 0         , 1 ];
        A12T = [cos(phi2_T) , -sin(phi2_T) , -L*sin(phi2_T); ...
                sin(phi2_T) ,  cos(phi2_T) , L*cos(phi2_T); ...
               0          , 0         , 1 ];
        A01T = [cos(phi1_T) , -sin(phi1_T) , -L*sin(phi1_T) ; ...
                sin(phi1_T) ,  cos(phi1_T) , L*cos(phi1_T); ...
                0        , 0          , 1  ];
        r4T = A34T;
        r4T = r4T(:,3);
        r4T = r4T(1:2);
        r3T = A34T*A23T;
        r3T = r3T(:,3);
        r3T = r3T(1:2);
        r2T = A34T*A23T*A12T;   %position link 2
        r2T = r2T(:,3);
        r2T = r2T(1:2);
        r1T = A34T*A23T*A12T*A01T;   %position link 1
        r1T = r1T(:,3);
        r1T = r1T(1:2);
    end
    R_all(:,i) = [r1T;r2T;r3T;r4T;state];
end
end


function animation(Rpos, movecoor, polar, lineObj)
dim = size(Rpos); 
for i = 1:dim(2)
    index = 2*(floor(i/20))+1; % 10*2
    if  Rpos(9,i) == 114 % rear
        refx = movecoor(1,index);
        [F,Max] = rotmot(0.0001, polar(:,i), polar(:,i+1), 0.14, 0.04, 0.1, 9.81, 1);
        pause(Max);
        set(lineObj.h0T,'xdata',[refx refx+Rpos(1,i)+0.2],'ydata',[0 0])
        set(lineObj.h1T,'xdata',[refx refx+Rpos(1,i)],'ydata',[0 Rpos(2,i)])
        set(lineObj.h2T,'xdata',[refx+Rpos(1,i) refx+Rpos(3,i)],'ydata',[Rpos(2,i) Rpos(4,i)])
        set(lineObj.h3T,'xdata',[refx+Rpos(3,i) refx+Rpos(5,i)],'ydata',[Rpos(4,i) Rpos(6,i)])
        set(lineObj.h4T,'xdata',[refx+Rpos(5,i) refx+Rpos(7,i)],'ydata',[Rpos(6,i) Rpos(8,i)])
        set(lineObj.h5T,'xdata',[refx+Rpos(5,i) refx+Rpos(7,i)-0.2],'ydata',[Rpos(8,i) Rpos(8,i)])
        drawnow;
    elseif Rpos(9,i) == 102 % front
        refx = movecoor(2,index);
        [F,Max] = rotmot(0.0001, polar(:,i), polar(:,i+1), 0.14, 0.04, 0.1, 9.81, 1);
        pause(Max);
        set(lineObj.h5T,'xdata',[refx refx+Rpos(7,i)-0.2],'ydata',[0 0])
        set(lineObj.h4T,'xdata',[refx refx+Rpos(7,i)],'ydata',[0 Rpos(8,i)])
        set(lineObj.h3T,'xdata',[refx+Rpos(7,i) refx+Rpos(5,i)],'ydata',[Rpos(8,i) Rpos(6,i)])
        set(lineObj.h2T,'xdata',[refx+Rpos(5,i) refx+Rpos(3,i)],'ydata',[Rpos(6,i) Rpos(4,i)])
        set(lineObj.h1T,'xdata',[refx+Rpos(3,i) refx+Rpos(1,i)],'ydata',[Rpos(4,i) Rpos(2,i)])
        set(lineObj.h0T,'xdata',[refx+Rpos(3,i) refx+Rpos(1,i)+0.2],'ydata',[Rpos(2,i) Rpos(2,i)])
        drawnow;
    end
end

end

function [F,Max] = rotmot(tormax, ang_1, ang_2, L, s, M, g, T)
dim = size(ang_1); a = 0.3;
    % obtain angle
    phi11 = ang_1(1);phi21 = ang_1(2);phi31 = ang_1(3);
    phi12 = ang_2(1);phi22 = ang_2(2);phi32 = ang_2(3);
    t1 = M*g*L*(sin(phi11)+2*sin(phi11+phi21*0.5)*cos(phi21*0.5));
    I1 = M*(s^2/3+L^2+(2*L*cos(phi21/2))^2);
    t2 = M*g*L*(sin(phi11)+sin(phi31));
    I2 = M*(s^2/3+2*L^2);
    t3 = M*g*L*(sin(phi31)+2*sin(phi31+phi21/2)*cos(phi21/2));
    I3 = M*(s^2/3 + L^2 +(2*L*cos(phi21/2))^2);
    % angular acceleration
    alpha1 = (tormax-t1)/I1;
    alpha2 = (tormax-t2)/I2;
    alpha3 = (tormax-t3)/I3;
    % angular difference fraction
    d1 = phi12-phi11; f1 = d1/alpha1/a/(1-a);
    if f1< 0 
        f1 = 0;
    end
    d2 = phi22-phi21; f2 = d2/alpha2/a/(1-a);
    if f2< 0 
        f2 = 0;
    end
    d3 = phi32-phi31; f3 = d3/alpha3/a/(1-a);
    if f3< 0 
        f3 = 0;
    end
    F = [f1,f2,f3];
    Max = floor(max(F))/T;

end

function expandangle = expandang(ang)
n = 10; count = 1;
dim = size(ang);
leg = ang(4,:);
ang = ang(1:3,:);

expandangle = [];
expandleg = [];
for i = 1:dim(2)-1
    expandangle(:,count) = ang(:,i);
    count = count+1;
    diff = (ang(:,i+1) - ang(:,i))/n;
    for j = 1:n-1
        expandangle(:,count) = ang(:,i)+diff*j;
        count = count+1;
    end
end
expandangle(:,count) = ang(:,dim(2));
for k = 2:dim(2)
    expandleg((k-2)*10+1:(k-1)*10)= leg(k);
end
new_dim = size(expandleg);
expandleg(:,new_dim(2)+1) = leg(dim(2));

expandangle = [expandangle;expandleg];
end

