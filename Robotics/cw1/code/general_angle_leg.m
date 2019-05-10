clear all; close all; clc;
% function for general angle
x = 0.12; y = 0; L = 0.1;
phi = 100;
c21 = (x^2+(y-L)^2)*0.5/L^2 -0.5;  disp(c21);
c22 = cosd(phi)-sind(phi)*tand(phi);  disp(c22);
c23 = x*tand(phi) + (y-L)/L;  disp(c23);
disp(c21-c23/c22);
ph2 = acosd(c21-c23/c22);

c31 = sind(phi)*(y-L) -cosd(phi)*x;
c32 = sind(ph2)*(cosd(phi)*(y-L)+sind(phi)*x-L)/(1+cosd(ph2));
ph3 = asind(0.5*(c31-c32)/L);

ph1 = phi-ph2-ph3; 

x_p = L*(sind(ph1)+ sind(ph1+ph2)+ sind(ph1+ph2+ph3));
y_p = L*(1+cosd(ph1)+cosd(ph1+ph2)+cosd(ph1+ph2+ph3));

A10 = [1 , 0 ,  0; ...
   0 , 1 , L; ...
   0 , 0 , 1 ];

A21 = [cosd(ph1) , sind(ph1) , L*sind(ph1) ; ...
       -sind(ph1) ,  cosd(ph1) , L*cosd(ph1); ...
        0        , 0          , 1  ];

A32 = [cosd(ph2) , sind(ph2) , L*sind(ph2); ...
       -sind(ph2) ,  cosd(ph2) , L*cosd(ph2); ...
       0         , 0         , 1 ];
A43 = [cosd(ph3) , sind(ph3) , L*sind(ph3); ...
       -sind(ph3) ,  cosd(ph3) , L*cosd(ph3); ...
        0         , 0         , 1 ];   
A54 = [ 0,1,L;...
       -1,0,0;... 
        0,0,1];

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

sp1 = subplot(2,1,1);
hold on
axis equal
axis([-0.4 0.4 0 0.45])
%target line object
lineObj.m1 = line(0,0,'color','r','LineWidth',5,'Marker','o');
set(lineObj.m1,'xdata',x,'ydata',y)
%link line objects
lineObj.h1 = line(0,0,'color','k','LineWidth',2);
lineObj.h2 = line(0,0,'color','k','LineWidth',2);
lineObj.h3 = line(0,0,'color','k','LineWidth',2);
lineObj.h4 = line(0,0,'color','k','LineWidth',2);
lineObj.h5 = line(0,0,'color','k','LineWidth',2);
%joint line objects
lineObj.d1 = line(0,0,'color','k','LineWidth',5,'Marker','o');
lineObj.d2 = line(0,0,'color','k','LineWidth',5,'Marker','o');
lineObj.d3 = line(0,0,'color','k','LineWidth',5,'Marker','o');
lineObj.d4 = line(0,0,'color','k','LineWidth',5,'Marker','o');


%set target position and pause
set(lineObj.h1,'xdata',[0 r11(1)],'ydata',[0 r11(2)])
set(lineObj.h2,'xdata',[r11(1) r22(1)],'ydata',[r11(2) r22(2)])
set(lineObj.h3,'xdata',[r22(1) r33(1)],'ydata',[r22(2) r33(2)])
set(lineObj.h4,'xdata',[r33(1) r44(1)],'ydata',[r33(2) r44(2)])
set(lineObj.h5,'xdata',[r44(1) r55(1)],'ydata',[r44(2) r55(2)])
set(lineObj.d1,'xdata',r11(1),'ydata',r11(2))
set(lineObj.d2,'xdata',r22(1),'ydata',r22(2))
set(lineObj.d3,'xdata',r33(1),'ydata',r33(2))
hold off;