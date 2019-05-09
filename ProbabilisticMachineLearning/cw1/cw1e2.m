% cw1e2
% cw1 e
clear all; close all; clc;
load('cw1e.mat');
figure(1)
subplot(3,2,1)
mesh(reshape(x(:,1),11,11),reshape(x(:,2),11,11),reshape(y,11,11));
hold on 
xlabel('x1')
ylabel('x2')
zlabel('y')
hold off 
subplot(3,2,2)
mesh(reshape(x(:,1),11,11),reshape(x(:,2),11,11),reshape(y,11,11));
hold on 
xlabel('x1')
ylabel('x2')
zlabel('y')
hold off
colorbar;

xs = apxGrid('expand',{linspace(-3,3,100)',linspace(-3,3,100)'});

% Model1 --------------------------------------------------------- %
meanfunce1 = [];
likfunce1 = @likGauss;
covfunce1 = @covSEard;
hype1 = struct('mean',[],'cov',[0 0 0], 'lik', 0);
hyp2e1 = minimize(hype1, @gp, -100, @infGaussLik, meanfunce1,...
    covfunce1, likfunce1, x,y);
nlmle1 = gp(hyp2e1, @infGaussLik, meanfunce1,...
    covfunce1, likfunce1, x,y);
[nlze1 dnlze1] = gp(hyp2e1, @infGaussLik, meanfunce1,...
    covfunce1, likfunce1, x,y,xs);

fe1 = [nlze1+2*sqrt(dnlze1); flipdim(nlze1-2*sqrt(dnlze1),1)];


subplot(3,2,3)
mesh(reshape(xs(:,1),100,100),reshape(xs(:,2),100,100),...
    reshape(nlze1,100,100))
fe1reshape = reshape(fe1, 100, 100, 2);
hold on;
% mesh(reshape(xs(:,1),100,100),reshape(xs(:,2),100,100),...
%     fe1reshape(:,:,1))
% mesh(reshape(xs(:,1),100,100),reshape(xs(:,2),100,100),...
%     fe1reshape(:,:,2))
xlabel('x1')
ylabel('x2')
zlabel('predictive mean')
hold off;

subplot(3,2,4)
mesh(reshape(xs(:,1),100,100),reshape(xs(:,2),100,100),...
    reshape(nlze1,100,100))
hold on;
xlabel('x1')
ylabel('x2')
zlabel('predictive mean')
hold off;

% Model2 --------------------------------------------------------- %
meanfunce2 = [];
likfunce2 = @likGauss;
covfunce2 = {@covSum, {@covSEard,@covSEard}};
hype2 = struct('mean',[], 'cov',[0 0 0 0 0 0],'lik',0);
hype2.cov = 0.1*randn(6,1); % if not, hype2.cov = [cov1, cov2], 
% where cov1 = cov2
hyp2e2 = minimize(hype2, @gp, -100, @infGaussLik, meanfunce2,...
    covfunce2, likfunce2, x,y);
nlmle2 = gp(hyp2e2, @infGaussLik, meanfunce2,...
    covfunce2, likfunce2, x,y);
[nlz dnlz] = gp(hyp2e2, @infGaussLik, meanfunce2,...
    covfunce2, likfunce2, x,y,xs);

fe2 = [nlz+2*sqrt(dnlz); flipdim(nlz-2*sqrt(dnlz),1)];


fe2reshape = reshape(fe2,100,100,2);
subplot(3,2,5)
mesh(reshape(xs(:,1),100,100), reshape(xs(:,2),100,100),...
    reshape(nlz, 100,100))
hold on 
xlabel('x1')
ylabel('x2')
zlabel('predictive mean')
hold off 

% colorbar
% hold on;
% mesh(reshape(xs(:,1),100,100), reshape(xs(:,2),100,100),...
%     fe2reshape(:,:,1));
% mesh(reshape(xs(:,1),100,100), reshape(xs(:,2),100,100),...
%     fe2reshape(:,:,2));
% hold off;

subplot(3,2,6)
mesh(reshape(xs(:,1),100,100), reshape(xs(:,2),100,100),...
    reshape(nlz, 100,100))
hold on;
xlabel('x1')
ylabel('x2')
zlabel('predictive mean')
hold off;

% --------------------------------------------------------- %
figure(2)
subplot(3,1,1)
mesh(reshape(xs(:,1),100,100), reshape(xs(:,2),100,100),...
    reshape(nlz, 100,100)-reshape(nlze1,100,100))
hold on;
colorbar;
xlabel('x1')
ylabel('x2')
zlabel('predictive mean difference')
hold off;

subplot(3,1,2)
mesh(reshape(xs(:,1),100,100), reshape(xs(:,2),100,100),...
    fe2reshape(:,:,1)- fe1reshape(:,:,1));
hold on;
xlabel('x1')
ylabel('x2')
zlabel('predictive error difference')
hold off;

subplot(3,1,3)
mesh(reshape(xs(:,1),100,100), reshape(xs(:,2),100,100),...
    fe2reshape(:,:,2)- fe1reshape(:,:,2));
hold on;
xlabel('x1')
ylabel('x2')
zlabel('predictive error difference')
hold off;






