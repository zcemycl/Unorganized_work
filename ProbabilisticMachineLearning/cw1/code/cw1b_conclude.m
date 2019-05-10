% cw1b final conclusion
clear all; close all; clc;
load('cw1a.mat')
% test sample
xs = linspace(-4,4,1000)';

n = 2; m = 2; 

meanfunc = [];
covfunc = @covSEiso;
covfunc2 = @covSEiso;
ell = 0; sf = 0; p = 2;
sn = 0;
likfunc = @likGauss;

hyp = struct('mean',[], 'cov',[ell,sf],'lik',sn);
% nlml = gp(hyp, @infGaussLik, meanfunc,...
%     covfunc2, likfunc, x,y);
% [mu0 s20] = gp(hyp, @infGaussLik,meanfunc,...
%     covfunc2,likfunc, x,y,xs);

hyp2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc,...
    covfunc2, likfunc, x,y);
nlml2 = gp(hyp2, @infGaussLik, meanfunc,...
    covfunc2, likfunc, x,y);
[mu s2] = gp(hyp2, @infGaussLik,meanfunc,...
    covfunc2,likfunc, x,y,xs);

f = [mu+2*sqrt(s2); flipdim(mu-2*sqrt(s2),1)];
fill([xs; flipdim(xs,1)], f, [7 7 7]/8)
hold on 
plot(xs, mu);
plot(x, y, '+')
axis([-4 4 -2 3]);
legend({'post-train predictive error bar',...
        'post-train predictive mean',...
        'training set'},'location','bestoutside')
xlabel('input, x'); ylabel('output, y')
hold off

