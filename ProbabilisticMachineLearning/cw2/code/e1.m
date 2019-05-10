    % e1
clear all; close all; clc;
load tennis_data

Matches = size(G);
p_win_list = [];
for i = 1:107
    indices = find(G==i);
    indices_win = [];
    w_count = 0;
    dim = size(indices);
    for j = 1:dim(1)
        if indices(j) <= 1801
            w_count = w_count+1;
        elseif indices(j) > 1801
            w_count = w_count-1;
        end
    end
    p_win_list(i) = w_count/dim(1);
end

Ms = p_win_list';
cw2;



