% PartD
% clear all; clc;

% smallworld;
% gridworld;
cliffworld;
[v2, pi2, all2] = qLearning(model,1000,50000,...
                            'constant',0.1,...
                            'constant',0.1);
figure(2)
plotVP(v2,pi2,paramSet);