% part c
% clear all; clc;

% smallworld;
% gridworld;
cliffworld;
[v, pi, all] = sarsa(model,1000,50000,...
                    'constant',0.1,...
                    'constant',0.1);
figure(1)
plotVP(v,pi,paramSet);