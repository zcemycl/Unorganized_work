% decode someone code
%% joint initial position 
%20,110,20 
%% waking ... ... 
%35,95,35;PhiVec = [0, -pi/6, -pi/6, -pi/6]'; 
step_height = 0.5; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tx = [linspace(2.2,2.2,10),linspace(2.2,2,10)]; Ty = [linspace(0.2,0.7,10), linspace(0.7,0.2,10)]; 

JAnglesTR = []; 
for j  = 1:20
    tPos = [Tx(j),Ty(j)]';
    JAngles = pure_angle(Tx(j),Ty(j),'r');
end


% Front % assuming fixed which is not true  % reality gap 
Tx = [linspace(2.5,2.5,10),linspace(2.5,2.8,10)]; Ty = [linspace(0.2,0.7,10), linspace(0.7,0.2,10)]; 
JAnglesTF = [];
PhiVec = [0,deg2rad(30), deg2rad(30), deg2rad(30)]'; 
for j  = 1:20
    tPos = [Tx(j),Ty(j)]';
    JAngles = pure_angle(-Tx(j),Ty(j),'f');
end
