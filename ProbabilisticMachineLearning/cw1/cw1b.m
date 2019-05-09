clear all; close all; clc;
load('cw1a.mat')
% test sample
xs = linspace(-4,4,1000)';

n = 2; m = 2; 

meanfunc = [];
covfunc = @covSEiso;
covfunc2 = @covSEiso;
ell = -1; sf = 0; p = 2;
sn = 5;
likfunc = @likGauss;
nlmlist = []; ellist = [];
j = 1;
for ell = -10:1:10
% ell increase or decrease --> linearize
% ell - widen the bar, + narrow the bar (-1/0 are better)
hyp = struct('mean',[], 'cov',[ell,sf],'lik',log(sn));
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
nlmlist(j) = nlml2;
ellist(j) = ell;
j = j+1;
end

% % subplot(n,m,1)
% f = [mu+2*sqrt(s2); flipdim(mu-2*sqrt(s2),1)];
% fill([xs; flipdim(xs,1)], f, [7 7 7]/8)
% hold on 
% plot(xs, mu);
% plot(x, y, '+')
% axis([-4 4 -2 2]);
% legend({'post-train predictive error bar',...
%         'post-train predictive mean',...
%         'training set'},'location','bestoutside')
% xlabel('input, x'); ylabel('output, y')
% hold off

% Comparing different length-scale, noise and signal
% ellist = [];
% sflist = [];
% i = 1; j = 1;
% for sf = -50:1:10
%     for ell = -20:1:20
%         hypc = struct('mean',[], 'cov',[ell,sf],'lik',0);
%         hyp2c = minimize(hypc, @gp, -100, @infGaussLik, meanfunc,...
%             covfunc2, likfunc, x,y);
%         nlml2c = gp(hyp2c, @infGaussLik, meanfunc,...
%             covfunc2, likfunc, x,y);
%         nlmlist(j,i) = nlml2c; 
%         ellist(j,i) = ell;
%         sflist(j,i) = sf;
%         i = i+1;
%     end
%     j = j+1; 
% end

% subplot(n,m,2);
% mesh(ellist, sflist, nlmlist)
% subplot(n,m,3);
% plot(ellist(:,1), nlmlist(:,1))
% subplot(n,m,4);
% plot(sflist(:,1), nlmlist(:,1))
% ell = -1; % min. try different sf 5.5, 10, 100   5-5.5is critical
% sf = 0; % -50 flat -> 0 -> 5 -> 10 (max) -> 15 (error exists) sf<sn, to positive definite
% %(random error), close to 0, has no all positive definite
% snlist = [];
% for sn = 1:20 %0-3
%     hypc = struct('mean',[], 'cov',[ell,sf],'lik',log(sn));
%     hyp2c = minimize(hypc, @gp, -100, @infGaussLik, meanfunc,...
%     covfunc2, likfunc, x,y);
%     nlml2c = gp(hyp2c, @infGaussLik, meanfunc,...
%     covfunc2, likfunc, x,y);
% 
%     nlmlist(i) = nlml2c;
%     snlist(i) = sn;
%     
%     i = i+1;
%     
% end

% ell = -1; % min. try different sf 5.5, 10, 100   5-5.5is critical
% sf = 0; % -50 flat -> 0 -> 5 -> 10 (max) -> 15 (error exists) sf<sn, to positive definite
% sn = 10;
% %(random error), close to 0, has no all positive definite
% sflist = [];
% for sf = 1:10 %0-3
%     hypc = struct('mean',[], 'cov',[ell,sf],'lik',log(sn));
%     hyp2c = minimize(hypc, @gp, -100, @infGaussLik, meanfunc,...
%     covfunc2, likfunc, x,y);
%     nlml2c = gp(hyp2c, @infGaussLik, meanfunc,...
%     covfunc2, likfunc, x,y);
% 
%     nlmlist(i) = nlml2c;
%     sflist(i) = sf;
%     
%     i = i+1;
%     
% end
% 
% snlist = [];
% sflist = [];
% ell = 0;
% nlmlist = [];
% i = 1; j = 1;
% for sf = 0:10
%     for sn = 0:10
%         hypc = struct('mean',[], 'cov',[ell,sf],'lik',sn);
%         hyp2c = minimize(hypc, @gp, -100, @infGaussLik, meanfunc,...
%             covfunc2, likfunc, x,y);
%         nlml2c = gp(hyp2c, @infGaussLik, meanfunc,...
%             covfunc2, likfunc, x,y);
%         nlmlist(j,i) = nlml2c; 
%         snlist(j,i) = sn;
%         sflist(j,i) = sf;
%         i = i+1;
%     end
%     j = j+1; 
% end