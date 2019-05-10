%% Simulation 
% Locomotion problem 2
clear all; close all; clc;
%% Constant
% Robot info.
% structure
L_foot  = 0.1 ; L_leg   = 0.14;
front_x = 0   ; front_y = 0   ;
rear_x  = 0.1 ; rear_y  = 0   ;
max_sep = 0.28; min_sep = 0.1 ;
% motion
max_x   = -0.18; 

% Platform info
plat_sx = -0.5;
plat_ex = -1.5;
plat_h  = 0.05;
end_all = -1.05;

%% Q-learning
% possible states chosen
pos_x_mov = linspace(max_x, -max_x, 37);
pos_y_mov = linspace(0,2*plat_h,5);
[posstate,indices] = matrixform(pos_x_mov,pos_y_mov,...
    front_x,rear_x,front_x,rear_x);
[reffx1, refrx1] = summary(posstate,indices,front_x,rear_x);
[posstate2,indices2] = matrixform(pos_x_mov,pos_y_mov,...
    front_x,rear_x,reffx1,refrx1);
[reffx2, refrx2] = summary(posstate2,indices2,reffx1,refrx1);
[posstate3,indices3] = matrixform(pos_x_mov,pos_y_mov,...
    front_x,rear_x,reffx2,refrx2);
[reffx3, refrx3] = summary(posstate3,indices3,reffx2,refrx2);
[posstate4,indices4] = matrixform(pos_x_mov,pos_y_mov,...
    front_x,rear_x,reffx3,refrx3);
[reffx4, refrx4] = summary(posstate4,indices4,reffx3,refrx3);
% This shows it will get under the terrain
[posstate5,indices5] = matrixform(pos_x_mov,pos_y_mov,...
    front_x,rear_x,reffx4,refrx4);
[reffx5, refrx5] = summary(posstate5,indices5,reffx4,refrx4);
[posstate6,indices6] = matrixform(pos_x_mov,pos_y_mov,...
    front_x,rear_x,reffx5,refrx5);
[reffx6, refrx6] = summary(posstate6,indices6,reffx5,refrx5);
[posstate7,indices7] = matrixform(pos_x_mov,pos_y_mov,...
    front_x,rear_x,reffx6,refrx6);
[reffx7, refrx7] = summary(posstate7,indices7,reffx6,refrx6);
[posstate8,indices8] = matrixform(pos_x_mov,pos_y_mov,...
    front_x,rear_x,reffx7,refrx7);
[reffx8, refrx8] = summary(posstate8,indices8,reffx7,refrx7);
[posstate9,indices9] = matrixform(pos_x_mov,pos_y_mov,...
    front_x,rear_x,reffx8,refrx8);
[reffx9, refrx9] = summary(posstate9,indices9,reffx8,refrx8);
[posstate10,indices10] = matrixform(pos_x_mov,pos_y_mov,...
    front_x,rear_x,reffx9,refrx9);
[reffx10, refrx10] = summary(posstate10,indices10,reffx9,refrx9);
[posstate11,indices11] = matrixform(pos_x_mov,pos_y_mov,...
    front_x,rear_x,reffx10,refrx10);
[reffx11, refrx11] = summary(posstate11,indices11,reffx10,refrx10);
%% Animation
figure(1)
lineObj = animInit(plat_sx,plat_ex,plat_h);
hold on;
% steps illustration
% front
scatter([front_x, reffx1, reffx2, reffx3,...
    reffx4, reffx5, reffx6, reffx7, reffx8,reffx9,reffx10,reffx11], ...
    [zeros(1,5),0.05*ones(1,7)], 'r')
% rear
scatter([rear_x, refrx1, refrx2, refrx3,...
    refrx4, refrx5, refrx6, refrx7, refrx8,refrx9,refrx10,refrx11], ...
    [zeros(1,8),0.05*ones(1,4)], 'b')
% figure(2)
% lineObj = animInit(plat_sx,plat_ex,plat_h);
% hold on;
% % scatter(posstate(1:185,1)+front_x,posstate(1:185,2),'r');
% [rank3r,rank3f] = ranking(posstate,3);
% %scatter(rank3f(:,1)+front_x,rank3f(:,2),'r');
% set(lineObj.h1T,'xdata',[0 0],...
%     'ydata',[0 0.14])
% set(lineObj.h2T,'xdata',[0 0.08],...
%     'ydata',[0.14 0.2549])
% set(lineObj.h3T,'xdata',[0.08 0.16],...
%     'ydata',[0.2549 0.14])
% set(lineObj.h4T,'xdata',[0.16 0.16],...
%     'ydata',[0.14 0])
% 
% scatter(rank3r(:,1)+rear_x,rank3r(:,2),'b');
% %scatter(posstate(186:370,1)+rear_x,posstate(186:370,2),'b');
%% Functions
% Summary for choices
function [reffx, refrx] = summary(state,indices,frontx,rearx)
choices = state(indices,:);
index_miny = find(choices(:,2) == min(choices(:,2)));
reffx = frontx; refrx = rearx;
if choices(index_miny, 3) == 114
    refrx = refrx + choices(index_miny,1);
elseif choices(index_miny,3) == 102
    reffx = reffx + choices(index_miny,1);
end

end

function reward = reward_pos(form,reffx,refrx)
dim = size(form);
reward = zeros(dim(1),1);
possible = zeros(dim(1),1);
for i = 1:dim(1)
    reward(i,1) = -10*form(i,1);
    
    if form(i,3) == 114
        ext = sqrt((reffx-(form(i,1)+refrx))^2+(0-(form(i,2)+0))^2);
        if or(ext > 0.28, reffx>(form(i,1)+refrx))
            possible(i,1) = 0;
%             reward(i,1) = -10000;
        else
            possible(i,1) = 1;
        end
    elseif form(i,3) == 102
        ext = sqrt((form(i,1)+reffx-refrx)^2+(form(i,2)+0-(0))^2);
        if or(ext > 0.28, (reffx+form(i,1))>refrx)
            possible(i,1) = 0;
%             reward(i,1) = -10000;
        else
            possible(i,1) = 1;
        end
    end
end
reward = [form,reward,possible];
end

function reward = reward_loc(form,frontx,rearx,reffx,refrx)
% Add platform info later
x_plt = -0.5; y_plt = 0.05;
%-----------------------%
dim = size(form);
loc = zeros(dim(1),1);
for i = 1:dim(1)
    fx = frontx;
    rx = rearx;
    if form(i,3) == 114
        loc(i,1) = -1000*(form(i,1)+refrx-rx);
        if form(i,1)+refrx < x_plt
            if form(i,2)+0-0 < y_plt
                loc(i,1) = -10000;
            end
        end
    elseif form(i,3) == 102
        loc(i,1) = -1000*(form(i,1)+reffx-fx);
        if form(i,1)+reffx < x_plt
            if form(i,2)+0-0 < y_plt
                loc(i,1) = -10000;
            end
        end
    end 
end
reward = [form,loc];
end

function sumr = sum_reward(form)
dim = size(form);
su = zeros(dim(1),1);
for i = 1:dim(1)
    if form(i,5) == 1
        su(i,1) = form(i,4)+form(i,6);
    elseif form(i,5) == 0
        su(i,1) = -10000;
    end
end
sumr = [form,su];
end

function [form,indices] = matrixform(xstate,ystate,frontx,rearx,reffx,refrx)
dimx= size(xstate); dimy = size(ystate);
form = zeros(dimx(2)*dimy(2),2);
for i = 1:dimx(2)*dimy(2)
    if rem(i,5) == 0
        form(i-4:i,1) = xstate(i/5);
        form(i-4:i,2) = ystate';
    end
end
form = [form;form];
leg = zeros(2*dimx(2)*dimy(2),1);
leg(1:dimx(2)*dimy(2)) = 'r'; leg(dimx(2)*dimy(2)+1:2*dimx(2)*dimy(2)) = 'f';
form = [form, leg];
form = reward_pos(form,reffx,refrx);
form = reward_loc(form,frontx,rearx,reffx,refrx);
form = sum_reward(form);
% select the greatest reward
indices = find(form(:,7)==max(form(:,7)));
form(indices,:)
end


% DONE
function lineObj = animInit(plat_sx,plat_ex,plat_h)
hold on
axis equal
axis([-1.2 0.4 -0.05 0.4])
plot([plat_sx,plat_sx,plat_ex,plat_ex],...
        [0,plat_h,plat_h,0],...
         'k','LineWidth',3);
plot([-1, -1], [0, 0.4], 'k-', 'LineWidth', 3);

%link line objects
lineObj.h0T = line(0,0,'color','k','LineWidth',2);
lineObj.h1T = line(0,0,'color','k','LineWidth',2);
lineObj.h2T = line(0,0,'color','k','LineWidth',2);
lineObj.h3T = line(0,0,'color','k','LineWidth',2);
lineObj.h4T = line(0,0,'color','k','LineWidth',2);
lineObj.h5T = line(0,0,'color','k','LineWidth',2);
end 

function [rankr,rankf] = ranking(state,mode) 
max_val = max(state(:,mode));
indices = find(state(:,mode)==max_val);
dim = size(indices);
indices_r = [];indices_f = [];
j = 1; k = 1;
for i = 1:dim(1)
    if indices(i) >185
        indices_r(j) = indices(i);
        j = j+1;
    else
        indices_f(k) = indices(k);
        k = k+1;
    end
end
rankr = state(indices_r,:);
rankf = state(indices_f,:);
end

