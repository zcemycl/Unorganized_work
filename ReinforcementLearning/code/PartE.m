% Part E
PartC;
PartD;
% sarsa
% alpha = 0.1, ep = 0.5(more explore), 0.2
% ep = 0.05, 0.0001(more optimal)
% ep = 0.1, alpha = 0.01 (more float), 0.5 ,
% 1(can't converge, more fluctuation)

% ep = 1 too random

% qlearn
% ep = 0.1
% alpha = 1(more optimal), 0.5, 0.1(less stable)
% alpha = 0.1
% ep = 0.5 (converge earlier), 0.8 (more explore)

% 1st: alpha=ep=0.1 (best plot)
% 2nd: iter dependent (only qlearn > sarsa)
% 3rd: alpha=0.1, ep iter dependent
% 4th: alpha iter dependent, ep=0.1 (sarsa more fluctuation)

figure(3);
plot(all/100);
hold on; 
plot(all2/100);
legend('sarsa','qlearn');
ylim([-100,0]);
xlim([0,300]);
xlabel('\times 100 epochs');
ylabel('sum of rewards per eposide');