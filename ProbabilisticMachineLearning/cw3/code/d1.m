% d1
% Plot
clear all;close all; clc;
seed = 100;
randn('seed',seed);
bmm;

dim = size(docs_mixpro);
hold on 
for i = 1:dim(1)
    plot(docs_mixpro(i,:))
end
legend(num2str(linspace(1,20,20)'),'Location','bestoutside')
title(['random seed:',num2str(seed)])
xlabel('iteration')
ylabel('Probability')
xlim([1, 20])






