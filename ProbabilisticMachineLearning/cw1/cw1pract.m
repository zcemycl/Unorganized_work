% Probabilistic Machine Learning 
% Coursework 1 practices
clear all;
close all;
clc;

f = figure(3);
n = 2; m = 5;
load('cw1a.mat');
% ------------------------------------------------- %
% a
% ------------------------------------------------- %
meanfunc = [];
covfunc = @covSEiso;
likfunc = @likGauss;
% test sample
xs = linspace(-4,4,60)';

hyp = struct('mean',[], 'cov',[-0.5,0],'lik',0);
hyp2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc,...
    covfunc, likfunc, x,y);
[mu s2] = gp(hyp2, @infGaussLik,meanfunc,...
    covfunc,likfunc, x,y,xs);

f = [mu+2*sqrt(s2); flipdim(mu-2*sqrt(s2),1)];
subplot(n,m,1);

fill([xs; flipdim(xs,1)], f, [7 7 7]/8)
hold on 
plot(xs, mu); 
axis([-4 4 -2 2])
plot(x, y, '+');
hold off

% ------------------------------------------------- %
% b 
% ------------------------------------------------- %
% non-convex --> ell,sc
% sf^2 * exp(-(x-x')/(2*ell^2))
ellvec = [-2,-1,0,1,2];
scvec = [-2,-1,0,1,2];
% ell = -1; sc = -2;
subplot(n,m,2);
for ell = -2:2:0.1
    for sc = -2:2:0.1
        hypb = struct('mean',[],'cov',[ell,sc],'lik',0);
        hypb2 = minimize(hypb, @gp, -100, @infGaussLik, meanfunc,...
        covfunc, likfunc, x,y);
        [mub s2b] = gp(hypb2, @infGaussLik, meanfunc,...
            covfunc, likfunc, x,y,xs);
        fb = [mub+2*sqrt(s2b); flipdim(mub-2*sqrt(s2b),1)];
        fill([xs; flipdim(xs,1)], fb, [8 8 8]/8)
        
        hold on 
        plot(xs, mub); 
        plot(x, y, '+');
        axis([-4 4 -2 3])
        alpha(0.1)
    end
end

% ------------------------------------------------- %
% c
% ------------------------------------------------- %
covfunc2 = @covPeriodic;
ell = -1; sf = 0; p = 2;
likfunc2 = @likGauss;
meanfunc2 = [];
hypc = struct('mean',[], 'cov',[ell,p,sf],'lik',0);
hyp2c = minimize(hypc, @gp, -100, @infGaussLik, meanfunc2,...
    covfunc2, likfunc2, x,y);
[muc s2c] = gp(hyp2c, @infGaussLik,meanfunc2,...
    covfunc2,likfunc2, x,y,xs);

fc = [muc+2*sqrt(s2c); flipdim(muc-2*sqrt(s2c),1)];
subplot(n,m,3);

fill([xs; flipdim(xs,1)], fc, [7 7 7]/8)
hold on 
plot(xs, muc); 
axis([-4 4 -2 2])
plot(x, y, '+');
hold off

% ------------------------------------------------- %
% d
% ------------------------------------------------- %
xd = linspace(-5,5,200)';
covfuncd = {@covProd, {@covPeriodic,@covSEiso}};
meanfuncd = @meanLinear;
likfuncd = @likGauss; sn = 0.1; 
hypd = struct('mean',0.5,'cov',[-0.5;0;0;2;0],...
    'lik',log(sn));
% without 1e-6*eye(200), K is not positive definite;
K = feval(covfuncd{:},hypd.cov,xd) + 1e-6*eye(200);
mud = feval(meanfuncd,hypd.mean,xd);
yd = chol(K)'*gpml_randn(0.15, 200, 1) + mud + exp(hyp.lik)*gpml_randn(0.2, 200, 1);

subplot(n,m,4);
plot(xd, yd, '+')

% ------------------------------------------------- %
% e
% ------------------------------------------------- %
load('cw1e.mat');

% visualization
subplot(n,m,5);
mesh(reshape(x(:,1),11,11),reshape(x(:,2),11,11),...
    reshape(y,11,11));

% 1st model
% ------------------------------------------------- %
meanfunce1 = [];
likfunce1 = @likGauss;
covfunce1 = @covSEard;
hype1 = struct('mean',[], 'cov',[-1;0],'lik',0);
hyp2e1 = minimize(hype1, @gp, -100, @infGaussLik, meanfunce1,...
    covfunce1, likfunce1, x(:,2),y);
hyp2e1b = minimize(hyp2e1, @gp, -100, @infGaussLik, meanfunce1,...
    covfunce1, likfunce1, x(:,1),y);
[mue1 s2e1] = gp(hyp2e1b, @infGaussLik,meanfunce1,...
    covfunce1,likfunce1, x(:,2),y,xs);
[mue1b s2e1b] = gp(hyp2e1b, @infGaussLik,meanfunce1,...
    covfunce1,likfunce1, x(:,1),y,xs);
fe1 = [mue1+2*sqrt(s2e1); flipdim(mue1-2*sqrt(s2e1),1)];
fe1b = [mue1b+2*sqrt(s2e1b); flipdim(mue1b-2*sqrt(s2e1b),1)];
subplot(n,m,6);

fill([xs; flipdim(xs,1)], fe1, [7 7 7]/8)
hold on 
plot(xs, mue1); 
axis([-4,4,-3,4]);
plot(x(:,2), y, '+');
hold off

subplot(n,m,7);
fill([xs; flipdim(xs,1)], fe1b, [7 7 7]/8)
hold on 
plot(xs, mue1b); 
axis([-4,4,-3,4]);
plot(x(:,1), y, '+');
hold off

% 2nd model
% ------------------------------------------------- %
meanfunce2 = [];
likfunce2 = @likGauss;
covfunce2 = {@covSum, {@covSEard,@covSEard}};

xs = apxGrid('expand',{linspace(-3,3,100)',linspace(-3,3,100)'});
hype2 = struct('mean',[], 'cov',[0 0 0 0 0 0],'lik',0);
hype2.cov = 0.1*randn(6,1);
hyp2e2 = minimize(hype2, @gp, -100, @infGaussLik, meanfunce2,...
    covfunce2, likfunce2, x,y);

[nlz dnlz] = gp(hyp2e2, @infGaussLik, meanfunce2,...
    covfunce2, likfunce2, x,y,xs);
% [nlml] = gp(hyp2, @infGaussLik, meanfunce2,...
%     covfunce2, likfunce2, x,y);



fe2 = [nlz+2*sqrt(dnlz); flipdim(nlz-2*sqrt(dnlz),1)];
subplot(n,m,8);
% fill([xs(:,1); flipdim(xs,1)], fe2, [7 7 7]/8)
% hold on 
% mesh(reshape(xs(:,1),121,1),y);
% axis([-4,4,-3,4]);
% mesh(reshape(x(:,1),11,11),reshape(x(:,2),11,11),...
%     reshape(y,11,11));
hold off

% ------------------------------------------------- %
% % try 3d prediction
% [xs1,xs2] = meshgrid(xs,xs);
% subplot(n,m,8);
% meanfunce3 = [];
% likfunce3 = @likGauss;
% covfunce3 = {@covSum, {@covSEard,@covSEard}};
% hype3 = struct('mean',[], 'cov',0.1*randn(6,1),'lik',0);
% hyp2e3 = minimize(hype3, @gp, -100, @infGaussLik, meanfunce3,...
%     covfunce3, likfunce3, x,y);
% [mue3 s2e3] = gp(hyp2e3, @infGaussLik,meanfunce3,...
%     covfunce3,likfunce3, x,y,xs1);
% 
% fe3 = [mue3+2*sqrt(s2e3); flipdim(mue3-2*sqrt(s2e3),1)];
% 
% fill([xs; flipdim(xs,1)], fe3, [7 7 7]/8)
% hold on 
% plot(xs, mue3); 
% axis([-4,4,-3,4]);
% plot(x, y, '+');
% hold off


