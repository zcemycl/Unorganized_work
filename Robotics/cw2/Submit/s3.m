% s3
%% Simulation 
clear all; close all; clc;
%% Constant
% Robot info.
% structure
L_foot  = 0.1 ; L_leg   = 0.14;
front_x = 0   ; front_y = 0   ;
rear_x  = 0.1 ; rear_y  = 0   ;
% motion
max_x   = -0.08; 
% Platform info
plt = [0,0.05;-0.1,0.1;-0.2,0.15;...
    -0.3,0.2;-0.4,0.25;-0.5,0];

%% Auto - Calculation
% move coordinates 
movecoor1 = position(front_x,rear_x, plt, L_leg);
% Rpos = transformation(movecoor1(:,1), L_leg);
%% Animation
[lineObj,pltall] = animInit(plt);
animation(movecoor1, lineObj, L_leg);
disp(movecoor1(6:8,:)/pi*180)
%% Functions
% DONE
function [lineObj,plt] = animInit(plat)
hold on
axis equal
axis([-0.7 0.2 -0.05 0.6])
dim = size(plat);
plt = zeros(12,2);
for i = 1:dim(1)
    tmpx = plat(i,1);
    if i == 1
        tmp1 = 0;
    else
        tmp1 = plat(i-1,2);
    end
    tmp2 = plat(i,2);
    plt(2*i-1,:) = [tmpx,tmp1];
    plt(2*i,:) = [tmpx,tmp2];
end
plot(plt(:,1),plt(:,2),'k','LineWidth',3);

%link line objects
lineObj.h0T = line(0,0,'color','k','LineWidth',2);
lineObj.h1T = line(0,0,'color','k','LineWidth',2);
lineObj.h2T = line(0,0,'color','k','LineWidth',2);
lineObj.h3T = line(0,0,'color','k','LineWidth',2);
lineObj.h4T = line(0,0,'color','k','LineWidth',2);
lineObj.h5T = line(0,0,'color','k','LineWidth',2);
end 

function movecoor = position(frontx,rearx, plat, L_leg)
dim = size(plat);
movecoor = zeros(4,4*(dim(1)-1));
tmpfx = frontx;
tmprx = rearx;
tmpfy = 0;
tmpry = 0;
for i = 1:dim(1)-1

    xd = plat(i+1,1)-plat(i,1);
    yd = plat(i+1,2)-plat(i,2);
    
    if i == dim(1)-1
        yd = 0.05;
    end

    fourvec = verthori(xd,yd,tmpfx,tmprx,tmpfy,tmpry);
    tmpfx = fourvec(1,4);
    tmprx = fourvec(2,4);
    tmpfy = fourvec(3,4);
    tmpry = fourvec(4,4);
    
    movecoor(:,4*(i-1)+1:4*i) = fourvec;
end
movecoor = [[frontx,rearx,0,0]',movecoor];
feet = zeros(1,4*(dim(1)-1)+1);
movecoor = [movecoor;feet];
for i = 1:4*(dim(1)-1)+1
    if i == 1
        movecoor(5,i) = 'f';
    end
    if rem(i,4) == 0
        movecoor(5,i-2) = 'f';
        movecoor(5,i-1) = 'f';
        movecoor(5,i)   = 'r';
        movecoor(5,i+1) = 'r';
    end
end
% % Add all angles info
polar = pos_angle(movecoor, L_leg);
movecoor = [movecoor;polar];
% % ------------------------------------- %
end

function out = verthori(dx,dy,rfx,rrx,rfy,rry)
    out = zeros(4,4);
    out(1,1) = rfx;
    out(1,2:4) = rfx+dx;
    out(3,:) = rfy+dy;
    out(2,1:3)=rrx;
    out(2,4)=rrx+dx;
    out(4,1:2)=rry;
    out(4,3:4) = rry+dy;
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



