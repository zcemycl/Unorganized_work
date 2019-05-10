function manipulator_two_legs
clear all 
close all
clc

%% parameters
l1 = 0.1;
l2 = 0.1;
l3 = 0.1;
l4 = 0.1;

%% initial conditions
phi1 = 0;
phi2 = 0;
phi3 = 0; 

phiVec = [phi1;phi2;phi3];
lVec = [l1,l2,l3,l4];

%% initialise animation
lineObj = animInit;
animation([phiVec,phiVec], lVec, lineObj, [-1,1]);

%% Run main routine
for i = 1:1
    %pDes = [[0.15,0.05]';0;1]; %000
    %pDes = [[0.12,0]';0;1]; %-pi/2
    pDes = [ginput(1)';0;1];    %get user defined target              
    phiVecNew = JacInv(phiVec,lVec,pDes);   %inverse kinematics
    animation([phiVec,phiVecNew],lVec,lineObj,pDes) %animation

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Animation initialisation
function lineObj = animInit()
    
    figure(1)
    hold on
    axis equal
    axis([-0.4 0.4 0 0.4])

    %target line object
    lineObj.m1 = line(0,0,'color','r','LineWidth',5,'Marker','o');
    %link line objects
    lineObj.h1 = line(0,0,'color','k','LineWidth',2);
    lineObj.h2 = line(0,0,'color','k','LineWidth',2);
    lineObj.h3 = line(0,0,'color','k','LineWidth',2);
    lineObj.h4 = line(0,0,'color','k','LineWidth',2);

    

    %joint line objects
    lineObj.d1 = line(0,0,'color','k','LineWidth',5,'Marker','o');
    lineObj.d2 = line(0,0,'color','k','LineWidth',5,'Marker','o');
    lineObj.d3 = line(0,0,'color','k','LineWidth',5,'Marker','o');
    lineObj.d4 = line(0,0,'color','k','LineWidth',5,'Marker','o');

    

end
%% Animation
function animation(phiMat,lVec,lineObj,PDes)
    
    PhiVec = phiMat(:,2);
    PhiVecOld = phiMat(:,1);
    [A10,A21,A32,A43] = HomCoord(PhiVec,lVec); %calculate transformation matrices


    
    r11 = A10*[0;l1;0;1];   %position link 1
    r22 = A10*A21*[0;l2;0;1];   %position link 2
    r33 = A10*A21*A32*[0;l3;0;1];
    r44 = A10*A21*A32*A43*[0;l4;0;1];

    
    %set target point
    set(lineObj.m1,'xdata',PDes(1),'ydata',PDes(2))
    
    %number of animated frames from start to end position
    n = 1000;
 
    if PhiVecOld ~= PhiVec
        
        for k = 1:n
            %avoid multiple rotations
            PhiVec = mod(PhiVec,2*pi);  
            PhiVecOld = mod(PhiVecOld,2*pi);
            dPhi = PhiVec-PhiVecOld;
            dC = abs(dPhi)>pi;
            dPhi = dPhi - dPhi.*dC - sign(dPhi.*dC).*(2*pi*dC-abs(dPhi.*dC));

            %calculate intermediate positions
            phiVecT = PhiVecOld + dPhi*k/n;
            [A10T,A21T,A32T,A43T] = HomCoord(phiVecT,lVec);
            r1T = A10T*[0;l1;0;1];
            r2T = A10T*A21T*[0;l2;0;1];
            r3T = A10T*A21T*A32T*[0;l3;0;1];
            r4T = A10T*A21T*A32T*A43T*[0;l4;0;1];
   
            %set intermediate positions and pause
            set(lineObj.h1,'xdata',[0 r1T(1)],'ydata',[0 r1T(2)])
            set(lineObj.h2,'xdata',[r1T(1) r2T(1)],'ydata',[r1T(2) r2T(2)])
            set(lineObj.h3,'xdata',[r2T(1) r3T(1)],'ydata',[r2T(2) r3T(2)])
            set(lineObj.h4,'xdata',[r3T(1) r4T(1)],'ydata',[r3T(2) r4T(2)])
            
            set(lineObj.d1,'xdata',r1T(1),'ydata',r1T(2))
            set(lineObj.d2,'xdata',r2T(1),'ydata',r2T(2))
            set(lineObj.d3,'xdata',r3T(1),'ydata',r3T(2))
            
            pause(1/n)
            
            
        end
    end
    
    %set target position and pause
     set(lineObj.h1,'xdata',[0 r11(1)],'ydata',[0 r11(2)])
     set(lineObj.h2,'xdata',[r11(1) r22(1)],'ydata',[r11(2) r22(2)])
     set(lineObj.h3,'xdata',[r22(1) r33(1)],'ydata',[r22(2) r33(2)])
     set(lineObj.h4,'xdata',[r33(1) r44(1)],'ydata',[r33(2) r44(2)])
            
            
     set(lineObj.d1,'xdata',r11(1),'ydata',r11(2))
     set(lineObj.d2,'xdata',r22(1),'ydata',r22(2))
     set(lineObj.d3,'xdata',r33(1),'ydata',r33(2))
%     set(lineObj.d2,'xdata',r22(1),'ydata',r22(2))
    drawnow
end

%% Homogeneous coordinate transform
function [A10,A21,A32,A43] = HomCoord(phiVec,lVec)
    
    Phi1 = phiVec(1);
    Phi2 = phiVec(2);
    Phi3 = phiVec(3);

    L1 = lVec(1);
    L2 = lVec(2);
    L3 = lVec(3);
    L4 = lVec(4);

    A10 = [1 , 0 , 0 , 0; ...
           0 , 1 , 0 , 0; ...
           0 , 0 , 1 , 0;
           0 , 0 , 0 , 1];

    A21 = [cos(Phi1) , -sin(Phi1) , 0 , 0; ...
           sin(Phi1) ,  cos(Phi1) , 0 , L2; ...
           0         , 0         , 1 , 0;
           0         , 0         , 0 , 1];
       
    A32 = [cos(Phi2) , -sin(Phi2) , 0 , 0; ...
           sin(Phi2) ,  cos(Phi2) , 0 , L3; ...
           0         , 0         , 1 , 0;
           0         , 0         , 0 , 1];
    A43 = [cos(Phi3) , -sin(Phi3) , 0 , 0; ...
           sin(Phi3) ,  cos(Phi3) , 0 , L4; ...
           0         , 0         , 1 , 0;
           0         , 0         , 0 , 1];

end

%% Jacobian inverse method (inverse kinematics)
function PhiVec = JacInv(PhiVec,LVec,PDes)
    
L4 = LVec(4)+0.1;   %offset target
j = 0;
[A10,A21,A32,A43] = HomCoord(PhiVec,LVec);
rCurr = A10*A21*A32*A43*[0;L4;0;1]; %current end effector position


pErr = norm(PDes-rCurr);    %current position error
errThresh = 0.013;  %error threshold
 
%  PhiVec 
% LVec
    while pErr > errThresh
        J = Jac(PhiVec,LVec);      

        dx =  pinv(J)*(PDes(1:2)-rCurr(1:2)); %pinv Moore-Penrose pseudoinverse
        PhiVec = PhiVec + dx;

        [A10,A21,A32,A43] = HomCoord(PhiVec,LVec);

        rCurr = A10*A21*A32*A43*[0;L4;0;1];
        pErr = norm(PDes-rCurr);

        j = j+1;

        if j>500 %interrupt if position not below threshold after 1000 iterations

            %run optimisation routine to minimise position error
            lambda = fminsearch(@(x) JacErr(x,PhiVec,PDes,LVec),[0;0;0]);
            PhiVec = PhiVec + lambda;
            PhiVec
            break
        end

    end
end

%% Calculate Jacobian
function J = Jac(phiVec,lVec)
    
phi1 = phiVec(1);
phi2 = phiVec(2);
phi3 = phiVec(3);

l1 = lVec(1);
l2 = lVec(2);
l3 = lVec(3);
l4 = lVec(4);


% Jacobian
 J = [  +l2*cos(phi1) + l3*cos(phi1+phi2) + l4*cos(phi1+phi2+phi3),  +l3*cos(phi1+phi2)+ l4*cos(phi1+phi2+phi3), l4*cos(phi1+phi2+phi3); ...
       -l2*sin(phi1)-l3*sin(phi1+phi2)-l4*sin(phi1+phi2+phi3), -l3*sin(phi1+phi2)-l4*sin(phi1+phi2+phi3),-l4*sin(phi1+phi2+phi3)];   
end

%% Objective function error minimisation
function pErr = JacErr(lambda,phiVec,pDes,lVec)

l4 = lVec(4);
phi = phiVec + lambda;
[A10,A21,A32,A43] = HomCoord(phi,lVec);
rCurr = A10*A21*A32*A43*[0;l4;0;1];
pErr = norm(pDes-rCurr); 
end

end