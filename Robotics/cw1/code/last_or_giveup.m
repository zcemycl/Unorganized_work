% last
clear all; close all; clc;
x = 0.15; y = 0.05;
L = 0.1;

phi1 = 90;
phi2 = -30; %-23.5
phi3 = 60;
g = -9.81;%gravity
f = -0.0003;%friction
m = 0.4; %mass
% steps
n = 200;

lineObj = initialization;
dphi2 = 20/n;
for i = 1:n
    [r11,r22,r33,r44,r55] = legspos(L,phi1,...
                            phi2+i*dphi2,phi3);
    draw(r11,r22,r33,r44,r55,L,lineObj,0,0);
    if r55(2) <= 0
        break
    end
    
end

m = 10;
dphi3 =  10/m;
for i = 1:m
    [r11,r22,r33,r44,r55] = legspos(L,phi1,...
                            phi2+65*dphi2,phi3-i*dphi3);
    draw(r11,r22,r33,r44,r55,L,lineObj,0,0);
end

% torque to force 
[Fx, Fy] = torqueforce(0, 0, 1, 90,-23.5,60,0.1);
disp(Fx);disp(Fy);

% timestep
dt = 0.005;
ax = Fx/m; ay = Fy/m;
vx = 0; vy = 0;
sx = 0; sy = 0;
for i = 1:100
    % x-axis
    if and(i == 1 , vx == 0)
        ax = ax;
    elseif and(i > 1 , vx < 0)
        ax = ax - f;
    elseif and(i > 1 , vx >= 0)
        disp(i)
        break
    end
    vx = vx + ax*dt;
    sx = sx + vx*dt;
%     % y-axis
%     if and(i == 1 , r55(2) == 0)
%         ay = ay;
%         vy = vy + ay*dt;
%     elseif and(i > 1 , r55(2) <= 0)
%         ay = 0;
%         vy = 0;
%     else
%         ay = ay+g;
%         vy = vy + ay*dt;
%     end
%     
%     sy = sy + vy*dt;
    draw(r11,r22,r33,r44,r55,L,lineObj,sx,0)
end

function lineObj = initialization()
    figure(1)
    hold on
    axis equal
    axis([-1.2 0.4 0 0.45])
    plot([-1.2 -0.5 -0.5],[0.0005 0.0005 0],...
         'k','LineWidth',3);

    %link line objects
    lineObj.h0 = line(0,0,'color','k','LineWidth',2);
    lineObj.h1 = line(0,0,'color','k','LineWidth',2);
    lineObj.h2 = line(0,0,'color','k','LineWidth',2);
    lineObj.h3 = line(0,0,'color','k','LineWidth',2);
    lineObj.h4 = line(0,0,'color','k','LineWidth',2);
    lineObj.h5 = line(0,0,'color','k','LineWidth',2);
    %joint line objects
    lineObj.d1 = line(0,0,'color','k','LineWidth',5,'Marker','o');
    lineObj.d2 = line(0,0,'color','k','LineWidth',5,'Marker','o');
    lineObj.d3 = line(0,0,'color','k','LineWidth',5,'Marker','o');
end

function [r11,r22,r33,r44,r55] = legspos(L,phi1,phi2,phi3)
    A10 = [1 , 0 ,  0; ...
           0 , 1 , L; ...
           0 , 0 , 1 ];

    A21 = [cosd(phi1) , sind(phi1) , L*sind(phi1) ; ...
           -sind(phi1) ,  cosd(phi1) , L*cosd(phi1); ...
           0         , 0         , 1  ];

    A32 = [cosd(phi2) , +sind(phi2) , L*sind(phi2); ...
           -sind(phi2) ,  cosd(phi2) , L*cosd(phi2); ...
           0         , 0         , 1 ];
    A43 = [cosd(phi3) , sind(phi3) , L*sind(phi3); ...
           -sind(phi3) ,  cosd(phi3) , L*cosd(phi3); ...
           0         , 0         , 1 ];
    A54 = [0 , 1, L;
           -1, 0, 0;
            0, 0, 1];

    r11 = A10(:,3);   %position link 1
    r11 = r11(1:2);
    r22 = A10*A21;   %position link 2
    r22 = r22(:,3);
    r22 = r22(1:2);
    r33 = A10*A21*A32;
    r33 = r33(:,3);
    r33 = r33(1:2);
    r44 = A10*A21*A32*A43;
    r44 = r44(:,3);
    r44 = r44(1:2);
    r55 = A10*A21*A32*A43*A54;
    r55 = r55(:,3);
    r55 = r55(1:2);
end

function draw(r11,r22,r33,r44,r55,L,lineObj,dx,dy)
    set(lineObj.h0,'xdata',[0+dx L+dx],'ydata',[dy dy])
    set(lineObj.h1,'xdata',[0+dx r11(1)+dx],'ydata',[dy r11(2)+dy])
    set(lineObj.h2,'xdata',[r11(1)+dx r22(1)+dx],'ydata',[r11(2)+dy r22(2)+dy])
    set(lineObj.h3,'xdata',[r22(1)+dx r33(1)+dx],'ydata',[r22(2)+dy r33(2)+dy])
    set(lineObj.h4,'xdata',[r33(1)+dx r44(1)+dx],'ydata',[r33(2)+dy r44(2)+dy])
    set(lineObj.h5,'xdata',[r44(1)+dx r55(1)+dx],'ydata',[r44(2)+dy r55(2)+dy])       

    set(lineObj.d1,'xdata',r11(1)+dx,'ydata',r11(2)+dy)
    set(lineObj.d2,'xdata',r22(1)+dx,'ydata',r22(2)+dy)
    set(lineObj.d3,'xdata',r33(1)+dx,'ydata',r33(2)+dy)
    drawnow;
end

function [Fx, Fy] = torqueforce(t1, t2, t3, phi1,phi2,phi3,L)
    J = L*[cosd(phi1)+cosd(phi1+phi2)+cosd(phi1+phi2+phi3),...
            cosd(phi1+phi2)+cosd(phi1+phi2+phi3),...
            cosd(phi1+phi2+phi3);...
            -sind(phi1)-sind(phi1+phi2)-sind(phi1+phi2+phi3),...
            -sind(phi1+phi2)-sind(phi1+phi2+phi3),...
            -sind(phi1+phi2+phi3)];
    F = pinv(J')*[t1,t2,t3]';
    Fx = F(1); Fy = F(2);
end