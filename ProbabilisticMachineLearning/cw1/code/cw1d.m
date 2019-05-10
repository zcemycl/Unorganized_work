% cw1d
clear all; close all; clc;
load('cw1a.mat')
% test sample
xs = linspace(-4,4,1000)';

xd = linspace(-5,5,200)';
covfuncd = {@covProd, {@covPeriodic,@covSEiso}};
meanfuncd = @meanLinear;
likfuncd = @likGauss; sn = 0.1; 
% [period/flat-period/linear;period.in;a.in;a.in;a.increase]

% without 1e-6*eye(200), K is not positive definite;
% 1e-6 --> the noise goes up
magnitudelist = [1e-6, 1e-3, 0.1, 1];
magnitude = magnitudelist(1);
varlist = [0, 1, 5];
for i = 1:3
hypd = struct('mean',0.5,'cov',[-0.5;0;0;2;varlist(i)],...
    'lik',log(sn));
K = feval(covfuncd{:},hypd.cov,xd) + magnitude*eye(200);
mud = feval(meanfuncd,hypd.mean,xd);
yd = chol(K)'*gpml_randn(0.15, 200, 1) + mud + exp(hypd.lik)*gpml_randn(0.2, 200, 1);
hypd2 = minimize(hypd, @gp, -100, @infGaussLik, meanfuncd,...
    covfuncd, likfuncd, xd,yd);
nlml2d = gp(hypd2, @infGaussLik, meanfuncd,...
            covfuncd, likfuncd, x,y);

[mud s2d] = gp(hypd2, @infGaussLik,meanfuncd,...
    covfuncd,likfuncd, xd,yd,xs);
fd = [mud+2*sqrt(s2d); flipdim(mud-2*sqrt(s2d),1)];
subplot(1,3,i)
fill([xs; flipdim(xs,1)], fd, [7 7 7]/8)
hold on 
plot(xs, mud);
plot(xd, yd, '+')
title(['std_f2:', num2str(varlist(i))])
% axis([-4 4 -6 6]);
legend({'post-train predictive error bar',...
        'post-train predictive mean',...
        'training set'},'location','south')
xlabel('input, x'); ylabel('output, y')
xlim([-4 4])
hold off

end
% ---------------------------------------------------------- %
% d_new1
% for i = 1:4
% 
% magnitude = magnitudelist(i);
% K = feval(covfuncd{:},hypd.cov,xd) + magnitude*eye(200);
% mud = feval(meanfuncd,hypd.mean,xd);
% yd = chol(K)'*gpml_randn(0.15, 200, 1) + mud + exp(hypd.lik)*gpml_randn(0.2, 200, 1);
% hypd2 = minimize(hypd, @gp, -100, @infGaussLik, meanfuncd,...
%     covfuncd, likfuncd, xd,yd);
% nlml2d = gp(hypd2, @infGaussLik, meanfuncd,...
%             covfuncd, likfuncd, x,y);
% 
% [mud s2d] = gp(hypd2, @infGaussLik,meanfuncd,...
%     covfuncd,likfuncd, xd,yd,xs);
% fd = [mud+2*sqrt(s2d); flipdim(mud-2*sqrt(s2d),1)];
% subplot(2,2,i)
% fill([xs; flipdim(xs,1)], fd, [7 7 7]/8)
% hold on 
% plot(xs, mud);
% plot(xd, yd, '+')
% title(['Magnitude of diagonal matrix:', num2str(magnitudelist(i))])
% axis([-4 4 -6 6]);
% legend({'post-train predictive error bar',...
%         'post-train predictive mean',...
%         'training set'},'location','south')
% xlabel('input, x'); ylabel('output, y')
% hold off
% 
% end
% ---------------------------------------------------------------- %
% d_new 2
% nlml2dlist = []; magnitudelist = [];
% j = 1;
% for magnitude = 1e-6:1e-2:1
%     K = feval(covfuncd{:},hypd.cov,xd) + magnitude*eye(200);
%     mud = feval(meanfuncd,hypd.mean,xd);
%     yd = chol(K)'*gpml_randn(0.15, 200, 1) + mud + exp(hypd.lik)*gpml_randn(0.2, 200, 1);
%     hypd2 = minimize(hypd, @gp, -100, @infGaussLik, meanfuncd,...
%                      covfuncd, likfuncd, xd,yd);
%     nlml2d = gp(hypd2, @infGaussLik, meanfuncd,...
%                 covfuncd, likfuncd, x,y);
%     nlml2dlist(j) = nlml2d;
%     magnitudelist(j) = magnitude;
%     j = j+1;
% end




