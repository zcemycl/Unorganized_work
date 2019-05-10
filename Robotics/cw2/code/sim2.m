%% Simulation 
% Locomotion problem 1
clear all; close all; clc;
%% Constant
% Time frame
T = 0.1;
% gravity
g = 10;
% Robot info.
% structure
L_foot  = 0.1 ; L_leg   = 0.14;
front_x = 0   ; front_y = 0   ;
rear_x  = 0.16; rear_y  = 0   ;
max_sep = 0.28; min_sep = 0.1 ;
% motion
max_x   = -0.08; 


% Platform info
plat_sx = -0.3;
plat_stx = -0.4;
plat_ex = -1;
plat_h  = 0.05;
plat_h2 = 0.1;

%% Auto - Calculation
% Total distance before the platform
[idealsize_bp, idealno_bp]=no_steps(front_x,plat_sx,max_x);
% move coordinates 
movecoor1 = position(front_x,rear_x, plat_h, ...
    idealsize_bp, idealno_bp, L_leg);
%% Animation
lineObj = animInit(plat_sx,plat_stx,plat_ex,plat_h,plat_h2);
animation(movecoor1, lineObj, L_leg);
%% Functions

% DONE
function lineObj = animInit(plat_sx,plat_stx,plat_ex,plat_h,plat_h2)
hold on
axis equal
axis([-2 1 -0.05 0.35])
plot([plat_sx,plat_sx,plat_stx,plat_stx,plat_ex,plat_ex],...
        [0,plat_h,plat_h,plat_h2,plat_h2,0],...
         'k','LineWidth',3);

%link line objects
lineObj.h0T = line(0,0,'color','k','LineWidth',2);
lineObj.h1T = line(0,0,'color','k','LineWidth',2);
lineObj.h2T = line(0,0,'color','k','LineWidth',2);
lineObj.h3T = line(0,0,'color','k','LineWidth',2);
lineObj.h4T = line(0,0,'color','k','LineWidth',2);
lineObj.h5T = line(0,0,'color','k','LineWidth',2);
end 

% DONE
function [idealsize, idealno]= no_steps(start_x,end_x,max_step)
% All possible step distance
step_size = [];
step_no = [];
d = end_x-start_x;
count = 1;
for i = 0:-0.01:max_step
    step_size(count) = i;
    step_no(count) = d/i;
    count = count+1;
end
% Choose the step size with min no.
dim = size(step_no);
step_no_integer = [];
count2 = 1;
for j = 1:dim(2)
    if rem(step_no(j),1) == 0
        step_no_integer(count2) = step_no(j);
        count2 = count2+1;
    end
end
idealno = min(abs(step_no_integer));
ind = find(step_no==idealno);
idealsize = step_size(ind);

end

% Flatsurface motion
function movecoor = position(frontx,rearx, height, distance, num, L_leg)
% ------------------------------------- %
% Generate front and rear legs position
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

% Add info about which feet is moving
feet = [];
feet(1) = 'r';
for i = 2:dim(2)
    diff = movecoor(1:2,i)-movecoor(1:2,i-1);
    if diff(1) ~= 0
        feet(i) = 'f';
    elseif diff(2) ~= 0
        feet(i) = 'r';
    end
end
movecoor = [movecoor;feet];
% Add all angles info
polar = pos_angle(movecoor, L_leg);
movecoor = [movecoor;polar];
end

% DONE
function polar = pos_angle(coor, L)
x = coor(2,:)-coor(1,:);
y = coor(4,:)-coor(3,:);
phi2 = acos((x.^2 + y.^2)./2./L^2 - 1);
phi3 = asin(0.5.*(x./L + sin(phi2).*((y-L)./L +1)./(cos(phi2)+1)));
phi1 = pi - phi2 - phi3; 
polar = [phi1;phi2;phi3];
end

function animation(movecoor, lineObj, L)
dim = size(movecoor); 
for i = 1:dim(2)
        n = 50;
        
        if i == dim(2)
            diff8 = zeros(8,1);
        else
            diff = (movecoor(6:8,i+1) - movecoor(6:8,i))/n;
            diff8 = [zeros(5,1); diff];
        end
        for j = 1:n-40
            tmp = movecoor(:,i)+diff8*j;
            Rpos = transformation(tmp, L);
            if  movecoor(5,i) == 114 % rear
                refx = movecoor(1,i);
                set(lineObj.h0T,'xdata',[refx refx+Rpos(1)+0.1],...
                    'ydata',[0 0])
                set(lineObj.h1T,'xdata',[refx refx+Rpos(1)],...
                    'ydata',[0 Rpos(2)])
                set(lineObj.h2T,'xdata',[refx+Rpos(1) refx+Rpos(3)],...
                    'ydata',[Rpos(2) Rpos(4)])
                set(lineObj.h3T,'xdata',[refx+Rpos(3) refx+Rpos(5)],...
                    'ydata',[Rpos(4) Rpos(6)])
                set(lineObj.h4T,'xdata',[refx+Rpos(5) refx+Rpos(7)],...
                    'ydata',[Rpos(6) Rpos(8)])
                set(lineObj.h5T,'xdata',[refx+Rpos(5) refx+Rpos(7)-0.1],...
                    'ydata',[Rpos(8) Rpos(8)])
                drawnow;
            elseif movecoor(5,i) == 102 % front
                refx = movecoor(2,i);
                set(lineObj.h5T,'xdata',[refx refx+Rpos(7)-0.1],...
                    'ydata',[0 0])
                set(lineObj.h4T,'xdata',[refx refx+Rpos(7)],...
                    'ydata',[0 Rpos(8)])
                set(lineObj.h3T,'xdata',[refx+Rpos(7) refx+Rpos(5)],...
                    'ydata',[Rpos(8) Rpos(6)])
                set(lineObj.h2T,'xdata',[refx+Rpos(5) refx+Rpos(3)],...
                    'ydata',[Rpos(6) Rpos(4)])
                set(lineObj.h1T,'xdata',[refx+Rpos(3) refx+Rpos(1)],...
                    'ydata',[Rpos(4) Rpos(2)])
                set(lineObj.h0T,'xdata',[refx+Rpos(3) refx+Rpos(1)+0.1],...
                    'ydata',[Rpos(2) Rpos(2)])
                drawnow;
            end
        end
end
end


function Rpos = transformation(movecoor, L)
phi1_T = movecoor(6); phi2_T = movecoor(7); 
phi3_T = movecoor(8); state = movecoor(5);
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
    Rpos = [r1T;r2T;r3T;r4T];
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
    Rpos = [r1T;r2T;r3T;r4T];
end

end