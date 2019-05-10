% Part b
clear all; clc;

% smallworld;
% gridworld;
cliffworld;
[v, pi] = policyIteration2(model,10000);
plotVP(v,pi,paramSet);