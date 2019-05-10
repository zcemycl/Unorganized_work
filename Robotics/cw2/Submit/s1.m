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
% after platform
front_x2 = front_x + idealsize_bp*idealno_bp - 0.1*2;
rear_x2 = front_x2 + min_sep;
[idealsize_bp2, idealno_bp2]=no_steps(front_x2,plat_ex,max_x);
% move coordinates 
movecoor1 = position(front_x,rear_x, plat_h, ...
    idealsize_bp, idealno_bp, L_leg, ...
    front_x2, rear_x2, plat_h2, idealsize_bp2, idealno_bp2);
Rpos = transformation(movecoor1(:,1), L_leg);
%% Animation
lineObj = animInit(plat_sx,plat_stx,plat_ex,plat_h,plat_h2);
animation(movecoor1, lineObj, L_leg);
%% Functions

% DONE
function lineObj = animInit(plat_sx,plat_stx,plat_ex,plat_h,plat_h2)
hold on
axis equal
axis([-1.2 0.4 -0.05 0.4])
plot([plat_sx,plat_sx,plat_stx,plat_stx,plat_ex,plat_ex],...
        [0,plat_h,plat_h,plat_h2,plat_h2,0],...
         'k','LineWidth',3);
plot([-0.7, -0.7], [0, 0.4], 'k-', 'LineWidth', 3);

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
function movecoor = position(frontx,rearx, height, distance, num, L_leg, frontx2, rearx2, height2, dist2, num2)
% ------------------------------------- %
% Generate front and rear legs position Flat
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
% Steps 
movecoor2 = zeros(4,8);
feet2 = [];
for s = 1:8
    if s == 1
        movecoor2(:,s) = movecoor(:,dim(2));
    else
        movecoor2(:,s) = movecoor2(:,s-1);
    end
    if rem(s,2) == 0
        
        if rem(s,4) == 0
            feet2(s-1) = 'r';
            feet2(s)   = 'r';
            movecoor2(2,s) = movecoor2(2,s)-0.1; % 0.1 stair size
            movecoor2(3,s-3:s) = height*(s/4); 
            movecoor2(4,s-1:s) = height*(s/4);
        else 
            feet2(s-1) = 'f';
            feet2(s)   = 'f';
            movecoor2(1,s) = movecoor2(1,s)-0.1; % 0.1 stair size
            
        end
    end
end
movecoor2 = [movecoor2;feet2];
% Flat walk over the stair
movecoor3 = zeros(4,4*num2+3);
movecoor3(1,1) = frontx2;
movecoor3(2,1:3) = rearx2;
movecoor3(3:4,:) = height2;
dim3 = size(movecoor3);
for i = 2:dim3(2)
    if rem(i,2) == 0
        if rem(i,4) == 0
            % Rear %
            % arrange x
            movecoor3(2,i) = movecoor3(2,i-1)+0.5*dist2;
            movecoor3(2,i+1:i+3) = movecoor3(2,i-1)+dist2;
            % arrange y
            movecoor3(4,i) = movecoor3(4,i) + height;
        else
            % Front %
            % arrange x
            movecoor3(1,i) = movecoor3(1,i-1)+0.5*dist2;
            movecoor3(1,i+1:i+3) = movecoor3(1,i-1)+dist2;
            % arrange y
            movecoor3(3,i) = movecoor3(3,i) + height;
            
        end
    end
end
movecoor3 = movecoor3(:,1:dim(2));
% ------------------------------------- %
% Add info about which feet is moving
feet = [];
feet(1) = 'r';
feet3 = [];
feet3(1) = 'f';
for i = 2:dim(2)
    diff = movecoor(1:2,i)-movecoor(1:2,i-1);
    if diff(1) ~= 0
        feet(i) = 'f';
    elseif diff(2) ~= 0
        feet(i) = 'r';
    end
end
for p = 2:dim(2)
    diff = movecoor3(1:2,p)-movecoor3(1:2,p-1);
    if diff(1) ~= 0
        feet3(p) = 'f';
    elseif diff(2) ~= 0
        feet3(p) = 'r';
    end
end
movecoor3 = [movecoor3;feet3];
movecoor = [movecoor;feet];
movecoor = [movecoor,movecoor2,movecoor3];
% ------------------------------------- %
% Add all angles info
polar = pos_angle(movecoor, L_leg);
movecoor = [movecoor;polar];
% ------------------------------------- %
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
        n = 20;
        for j = 1:n
            tmp = movecoor(:,i);
            Rpos = transformation(tmp, L);
            if  movecoor(5,i) == 114 % rear
                refx = movecoor(1,i);
                refy = movecoor(3,i);
                set(lineObj.h0T,'xdata',[refx refx+Rpos(1)+0.1],...
                    'ydata',[refy refy])
                set(lineObj.h1T,'xdata',[refx refx+Rpos(1)],...
                    'ydata',[refy refy+Rpos(2)])
                set(lineObj.h2T,'xdata',[refx+Rpos(1) refx+Rpos(3)],...
                    'ydata',[refy+Rpos(2) refy+Rpos(4)])
                set(lineObj.h3T,'xdata',[refx+Rpos(3) refx+Rpos(5)],...
                    'ydata',[refy+Rpos(4) refy+Rpos(6)])
                set(lineObj.h4T,'xdata',[refx+Rpos(5) refx+Rpos(7)],...
                    'ydata',[refy+Rpos(6) refy+Rpos(8)])
                set(lineObj.h5T,'xdata',[refx+Rpos(5) refx+Rpos(7)-0.1],...
                    'ydata',[refy+Rpos(8) refy+Rpos(8)])
                drawnow;
            elseif movecoor(5,i) == 102 % front
                refx = movecoor(2,i);
                refy = movecoor(4,i);
                set(lineObj.h5T,'xdata',[refx refx+Rpos(7)-0.1],...
                    'ydata',[refy refy])
                set(lineObj.h4T,'xdata',[refx refx+Rpos(7)],...
                    'ydata',[refy refy+Rpos(8)])
                set(lineObj.h3T,'xdata',[refx+Rpos(7) refx+Rpos(5)],...
                    'ydata',[refy+Rpos(8) refy+Rpos(6)])
                set(lineObj.h2T,'xdata',[refx+Rpos(5) refx+Rpos(3)],...
                    'ydata',[refy+Rpos(6) refy+Rpos(4)])
                set(lineObj.h1T,'xdata',[refx+Rpos(3) refx+Rpos(1)],...
                    'ydata',[refy+Rpos(4) refy+Rpos(2)])
                set(lineObj.h0T,'xdata',[refx+Rpos(3) refx+Rpos(1)+0.1],...
                    'ydata',[refy+Rpos(2) refy+Rpos(2)])
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