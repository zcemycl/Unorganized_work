clear all; clc;
% part a

% smallworld;
% gridworld;
cliffworld;
[v, pi] = valueIteration(model,10000);
plotVP(v,pi,paramSet);

