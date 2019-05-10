clear all;

%% joint initial position 
%20,110,20 
%% waking ... ... 
%35,95,35;PhiVec = [0, -pi/6, -pi/6, -pi/6]'; 
step_height = 0.5; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tx = [linspace(2.2,2.2,10),linspace(2.2,2,10)]; Ty = [linspace(0.2,0.7,10), linspace(0.7,0.2,10)]; 
%Rear
JAnglesTR = []; 
for j  = 1:20
    tPos = [Tx(j),Ty(j),1]';
    [JAngles, PhiVec] = JacInvR(PhiVec,LVec,tPos,"rear");
    JAnglesTR = [JAnglesTR,JAngles];
    disp(j);
end
% Front % assuming fixed which is not true  % reality gap 
Tx = [linspace(2.5,2.5,10),linspace(2.5,2.8,10)]; Ty = [linspace(0.2,0.7,10), linspace(0.7,0.2,10)]; 
JAnglesTF = [];
PhiVec = [0,deg2rad(30), deg2rad(30), deg2rad(30)]'; 
for j  = 1:20
    tPos = [Tx(j),Ty(j),1]';
    [JAngles, PhiVec] = JacInvR(PhiVec,LVec,tPos,"front");
    JAnglesTF = [JAnglesTF,JAngles];
    disp(j);
end
%%%%%%%%% Walking ?? %%%%%%%%%
pause(2)
for i = 1:20 
for j = 1:20
    angles(JAnglesTR(2,j),JAnglesTR(3,j),JAnglesTR(4,j), cal);
end
pause(0.5)
for j = 1:19
    angles(JAnglesTF(4,j), JAnglesTF(3,j),JAnglesTF(2,j), cal);
end 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Walking + Obstacle + Walking 
th1 = [linspace(20,30,10), linspace(30,30,10), linspace(30,30,10), linspace(30,10,10), linspace(10,10,10), linspace(10,20,10)]; 
th2 = [linspace(110,100,10), linspace(100,95,10), linspace(95,95,10), linspace(95,95,10), linspace(95,100,10), linspace(100,110,10)]; 
th3 = [linspace(20,10,10), linspace(10,10,10), linspace(10,30,10), linspace(30,30,10), linspace(30,30,10), linspace(30,20,10)]; 
angles(20,115,20, cal);
pause(3)
for j = 1:15
for i = 1:60
    angles(th1(i),th2(i), th3(i), cal); 
end
end
% Obstacle 
th1 = [linspace(20,40,10), linspace(40,60,10),linspace(60,20,20),linspace(20,20,10),linspace(20,40,10),90, linspace(90,70,10), linspace(70,70,10),linspace(70,20,10)]; 
th2 = [linspace(110,110,10), linspace(110,60,10),linspace(60,60,20),linspace(60,60,10),linspace(60,60,10),0,linspace(0,0,10), linspace(0,40,10), linspace(40,110,10)];
th3 = [linspace(20,-10,10), linspace(-10,-10,10), linspace(-10,0,20),linspace(0,70,10),linspace(70,70,10),90, linspace(90,80,10), linspace(80,80,10), linspace(80,20,10)]; 
for i = 1:length(th1)
    angles(th1(i),th2(i),th3(i), cal)
    pause(0.05)
end

th1 = [linspace(20,30,10), linspace(30,30,10), linspace(30,30,10), linspace(30,10,10), linspace(10,10,10), linspace(10,20,10)]; 
th2 = [linspace(110,100,10), linspace(100,95,10), linspace(95,95,10), linspace(95,95,10), linspace(95,100,10), linspace(100,110,10)]; 
th3 = [linspace(20,10,10), linspace(10,10,10), linspace(10,30,10), linspace(30,30,10), linspace(30,30,10), linspace(30,20,10)]; 
angles(20,115,20, cal);

for j = 1:15
for i = 1:60
    angles(th1(i),th2(i), th3(i), cal); 
end
end
%% Stairs

% steps 15 

th1 = [linspace(20,30,10), linspace(30,30,10), linspace(30,30,10), ...
    linspace(30,10,10), linspace(10,10,10), linspace(10,20,10)]; 
th2 = [linspace(110,100,10), linspace(100,95,10), linspace(95,95,10), ...
    linspace(95,95,10), linspace(95,100,10), linspace(100,110,10)]; 
th3 = [linspace(20,10,10), linspace(10,10,10), linspace(10,30,10), ...
    linspace(30,30,10), linspace(30,30,10), linspace(30,20,10)]; 
angles(20,115,20, cal);
pause(2)
for j = 1:15
for i = 1:60
    angles(th1(i),th2(i), th3(i), cal); 
end
end

% stairs 

th1 = [linspace(20,40,10),linspace(40,60,10), linspace(60,60,10), linspace(60,55,10), ...
    linspace(55,50,10), linspace(50,50,10), linspace(50,35,10), linspace(35,35,10), ...
    linspace(35,35,10), linspace(35,30,10), linspace(30,30,10), linspace(30,20,10),...
    linspace(20,10,10), linspace(10,-10,10), linspace(-10,-10,10), linspace(-10,20,10),...
    linspace(20,30,10), linspace(30,60,10), linspace(60,60,10),linspace(50,50,5),linspace(45,45,5),linspace(45,45,10), ....
    linspace(45,60,10), linspace(60,50,10), linspace(50,35,10), linspace(35,35,10), ...
    linspace(35,35,10), linspace(35,30,10), linspace(30,30,10), linspace(30,20,10),...
    linspace(20,10,10), linspace(10,-10,10), linspace(-10,-10,10), linspace(-10,20,10)]; 
%
th2 = [linspace(110,110,10), linspace(110,75,10), linspace(75,60,10), linspace(60,55,10), ...
    linspace(55,65,10), linspace(65,73,10), linspace(73,85,10), linspace(85,95,10),...
    linspace(95,100,10), linspace(100,110,10), linspace(110,120,10), linspace(120,120,10),...
    linspace(120,120,10), linspace(120,120,10), linspace(120,110,10),linspace(110,110,10),...
    linspace(110,110,10), linspace(110,100,10),linspace(100,70,10),linspace(90,90,5),linspace(100,100,5),linspace(100,80,10),...
    linspace(80,60,10), linspace(60,73,10), linspace(73,85,10), linspace(85,95,10),...
    linspace(95,100,10), linspace(100,110,10), linspace(110,120,10), linspace(120,120,10),...
    linspace(120,120,10), linspace(120,120,10), linspace(120,110,10),linspace(110,110,10)]; 
%
th3 = [linspace(20,-10,10), linspace(-10,-10,10), linspace(-10,40,10), linspace(40,30,10), ...
    linspace(30,20,10), linspace(20,10,10), linspace(10,-5,10), linspace(-5,-20,10),...
    linspace(-20,-30,10), linspace(-30,-45,10), linspace(-45,-45,10), linspace(-45,-35,10)...
    linspace(-35,-20,10),linspace(-20,0,10), linspace(0,20,10), linspace(20,20,10),...
    linspace(20,0,10),linspace(0,-10,10),linspace(-10,25,10),linspace(10,10,5),linspace(0,0,5),linspace(0,-20,10),...
    linspace(-20,35,10), linspace(35,10,10), linspace(10,-5,10), linspace(-5,-20,10),...
    linspace(-20,-30,10), linspace(-30,-45,10), linspace(-45,-45,10), linspace(-45,-35,10)...
    linspace(-35,-20,10),linspace(-20,0,10), linspace(0,20,10), linspace(20,20,10)]; 

for i = 1:length(th1)
    angles(th1(i),th2(i),th3(i), cal)
    pause(0.05)
end

% steps 4

th1 = [linspace(20,30,10), linspace(30,30,10), linspace(30,30,10), linspace(30,10,10), linspace(10,10,10), linspace(10,20,10)]; 
th2 = [linspace(110,100,10), linspace(100,95,10), linspace(95,95,10), linspace(95,95,10), linspace(95,100,10), linspace(100,110,10)]; 
th3 = [linspace(20,10,10), linspace(10,10,10), linspace(10,30,10), linspace(30,30,10), linspace(30,30,10), linspace(30,20,10)]; 
angles(20,115,20, cal);
pause(3)
for j = 1:4
for i = 1:60
    angles(th1(i),th2(i), th3(i), cal); 
end
end

