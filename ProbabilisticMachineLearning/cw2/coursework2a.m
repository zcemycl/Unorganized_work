% coursework part a
n = 4; m = 4;
% x_cov 
c_all = [];
for i = 1:107
[c, lg] = xcov(w_all(i,:), 100, 'coeff');
c_all(i,:) = c;
end

c_sum = [];
c_sum = sum(c_all, 2);

% Sample player 16
index_16_4 = [];
ind_16_4 = [];
for i = 1:275
    index_16_4(i) = 4*i;
    ind_16_4(i) = w_all(16,4*i);
end

% ---------------------------------------- %
figure(1)
title('Gibbs Sampler')
subplot(n,m,1)
plot(w_all(16,:))
hold on 
plot(w_sum_all(16,:)); 
title(['1st:', W(16)])
xlim([0,1100]);
subplot(n,m,2)
plot(var_all(16,:));
title(['Variance'])
xlim([0,1100]);
subplot(n,m,3)
plot(index_16_4, ind_16_4);
title(['Every 4th sample'])
xlim([0,1100]);

index_16_5 = [];
ind_16_5 = [];
for i = 1:220
    index_16_5(i) = 5*i;
    ind_16_5(i) = w_all(16,5*i);
end
subplot(n,m,4)
plot(index_16_5, ind_16_5);
title(['Every 5th sample'])
xlim([0,1100]);


% Sample player 1
index_1_2 = [];
ind_1_2 = [];
for i = 1:550
    index_1_2(i) = 2*i;
    ind_1_2(i) = w_all(1,2*i);
end
subplot(n,m,5)
plot(w_all(1,:))
hold on 
plot(w_sum_all(1,:)); 
title(['2nd:',W(1)])
xlim([0,1100]);
subplot(n,m,6)
title(['Variance'])
plot(var_all(1,:));
xlim([0,1100]);
subplot(n,m,7)
plot(index_1_2, ind_1_2);
title(['Every 2nd sample'])
xlim([0,1100]);
subplot(n,m,8)
plot(index,g_sample(1,:))
title(['Every 5th sample'])
xlim([0,1100]);

% Sample player 5
index_5_2 = [];
ind_5_2 = [];
for i = 1:220
    index_5_2(i) = 5*i;
    ind_5_2(i) = w_all(1,5*i);
end
subplot(n,m,9)
plot(w_all(5,:))
hold on 
plot(w_sum_all(5,:)); 
title(['3rd:',W(5)])
xlim([0,1100]);
subplot(n,m,10)
title(['Variance'])
plot(var_all(5,:));
xlim([0,1100]);
subplot(n,m,11)
plot(index_5_2, ind_5_2);
title(['Every 5th sample'])
xlim([0,1100]);
subplot(n,m,12)
plot(index_5_2, ind_5_2);
title(['Every 5th sample'])
xlim([0,1100]);

% Sample player 77
index_77_2 = [];
ind_77_2 = [];
for i = 1:1100
    index_77_2(i) = i;
    ind_77_2(i) = w_all(1,i);
end
subplot(n,m,13)
plot(w_all(77,:))
hold on 
plot(w_sum_all(77,:)); 
title(['53th:',W(77)])
xlim([0,1100]);
subplot(n,m,14)
title(['Variance'])
plot(var_all(77,:));
xlabel('iteration')
xlim([0,1100]);
subplot(n,m,15)
plot(index_77_2, ind_77_2);
title(['Every 1st sample'])
xlim([0,1100]);
subplot(n,m,16)
plot(index, g_sample(5,:));
title(['Every 5th sample'])
xlim([0,1100]);

% ---------------------------------------- %
figure(2)
plot(c_sum)
xlim([1,107])
ylabel('N-th sample kept')
xlabel('Player index')
