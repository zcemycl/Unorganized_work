% cw1 for each
clear all;
close all;
clc;

f = figure(3);
load('cw1a.mat');
% ------------------------------------------------- %
% Optimization
% ------------------------------------------------- %
meanfunc = [];
covfunc = @covSEiso;
covfunc2 = @covPeriodic;
ell = -1; sf = 0; p = 2;
likfunc = @likGauss;
% test sample
xs = linspace(-4,4,1000)';
% ell increase or decrease --> linearize
% ell - widen the bar, + narrow the bar (-1/0 are better)
hyp = struct('mean',[], 'cov',[ell,p,sf],'lik',0);
nlml = gp(hyp, @infGaussLik, meanfunc,...
    covfunc2, likfunc, x,y);
[mu0 s20] = gp(hyp, @infGaussLik,meanfunc,...
    covfunc2,likfunc, x,y,xs);

hyp2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc,...
    covfunc2, likfunc, x,y);
nlml2 = gp(hyp2, @infGaussLik, meanfunc,...
    covfunc2, likfunc, x,y);
[mu s2] = gp(hyp2, @infGaussLik,meanfunc,...
    covfunc2,likfunc, x,y,xs);
% --------------------------------------------------- %
% b
nlmlstor = []; j = 1;
elllist = [];
for i = -10:0.1:10
    hyp = struct('mean',[], 'cov',[i,0],'lik',0);
    hyp2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc,...
    covfunc, likfunc, x,y);
    nlml2 = gp(hyp2, @infGaussLik, meanfunc,...
    covfunc, likfunc, x,y);
    
    nlmlstor(j) = nlml2;
    elllist(j) = i;
    j = j+1;
end
figure(3);
plot(elllist, nlmlstor);
% --------------------------------------------------- %
f0 = [mu0+2*sqrt(s20); flipdim(mu0-2*sqrt(s20),1)];
f = [mu+2*sqrt(s2); flipdim(mu-2*sqrt(s2),1)];
% % ---------------------------------------------- %
% hyp3 = struct('mean',[], 'cov',[-10,0],'lik',0);
% nlml3 = gp(hyp3, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% [mu3 s23] = gp(hyp3, @infGaussLik,meanfunc,...
%     covfunc,likfunc, x,y,xs);
% 
% hyp4 = minimize(hyp3, @gp, -100, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% nlml4 = gp(hyp4, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% [mu4 s24] = gp(hyp4, @infGaussLik,meanfunc,...
%     covfunc,likfunc, x,y,xs);
% f3 = [mu3+2*sqrt(s23); flipdim(mu3-2*sqrt(s23),1)];
% f4 = [mu4+2*sqrt(s24); flipdim(mu4-2*sqrt(s24),1)];
% % ---------------------------------------------- %
% hyp5 = struct('mean',[], 'cov',[-10,0],'lik',0);
% nlml5 = gp(hyp5, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% [mu5 s25] = gp(hyp5, @infGaussLik,meanfunc,...
%     covfunc,likfunc, x,y,xs);
% 
% hyp6 = minimize(hyp5, @gp, -100, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% nlml6 = gp(hyp6, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% [mu6 s26] = gp(hyp6, @infGaussLik,meanfunc,...
%     covfunc,likfunc, x,y,xs);
% f5 = [mu5+2*sqrt(s25); flipdim(mu5-2*sqrt(s25),1)];
% f6 = [mu6+2*sqrt(s26); flipdim(mu6-2*sqrt(s26),1)];
% % ----------------------------------------------------- %
% hyp7 = struct('mean',[], 'cov',[-0.25,-5],'lik',0);
% nlml7 = gp(hyp7, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% [mu7 s27] = gp(hyp7, @infGaussLik,meanfunc,...
%     covfunc,likfunc, x,y,xs);
% 
% hyp28 = minimize(hyp7, @gp, -100, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% nlml28 = gp(hyp28, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% [mu8 s28] = gp(hyp28, @infGaussLik,meanfunc,...
%     covfunc,likfunc, x,y,xs);
% f7 = [mu7+2*sqrt(s27); flipdim(mu7-2*sqrt(s27),1)];
% f8 = [mu8+2*sqrt(s28); flipdim(mu8-2*sqrt(s28),1)];
% % ---------------------------------------------- %
% hyp9 = struct('mean',[], 'cov',[-10,5],'lik',0);
% nlml9 = gp(hyp9, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% [mu9 s29] = gp(hyp9, @infGaussLik,meanfunc,...
%     covfunc,likfunc, x,y,xs);
% 
% hyp10 = minimize(hyp9, @gp, -100, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% nlml10 = gp(hyp10, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% [mu10 s210] = gp(hyp10, @infGaussLik,meanfunc,...
%     covfunc,likfunc, x,y,xs);
% f9 = [mu9+2*sqrt(s29); flipdim(mu9-2*sqrt(s29),1)];
% f10 = [mu10+2*sqrt(s210); flipdim(mu10-2*sqrt(s210),1)];
% % ---------------------------------------------- %
% hyp11 = struct('mean',[], 'cov',[-10,-5],'lik',0);
% nlml11 = gp(hyp11, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% [mu11 s211] = gp(hyp11, @infGaussLik,meanfunc,...
%     covfunc,likfunc, x,y,xs);
% 
% hyp12 = minimize(hyp11, @gp, -100, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% nlml12 = gp(hyp12, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% [mu12 s212] = gp(hyp12, @infGaussLik,meanfunc,...
%     covfunc,likfunc, x,y,xs);
% f11 = [mu11+2*sqrt(s25); flipdim(mu11-2*sqrt(s25),1)];
% f12 = [mu12+2*sqrt(s26); flipdim(mu12-2*sqrt(s26),1)];
% 
% % ---------------------------------------------- %
% hyp13 = struct('mean',[], 'cov',[-0.25,0],'lik',5);
% nlml13 = gp(hyp13, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% [mu13 s213] = gp(hyp13, @infGaussLik,meanfunc,...
%     covfunc,likfunc, x,y,xs);
% 
% hyp214 = minimize(hyp13, @gp, -100, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% nlml214 = gp(hyp214, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% [mu14 s214] = gp(hyp214, @infGaussLik,meanfunc,...
%     covfunc,likfunc, x,y,xs);
% f13 = [mu13+2*sqrt(s213); flipdim(mu13-2*sqrt(s213),1)];
% f14 = [mu14+2*sqrt(s214); flipdim(mu14-2*sqrt(s214),1)];
% % ---------------------------------------------- %
% hyp15 = struct('mean',[], 'cov',[-10,0],'lik',-5);
% nlml15 = gp(hyp15, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% [mu15 s215] = gp(hyp15, @infGaussLik,meanfunc,...
%     covfunc,likfunc, x,y,xs);
% 
% hyp16 = minimize(hyp15, @gp, -100, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% nlml16 = gp(hyp16, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% [mu16 s216] = gp(hyp16, @infGaussLik,meanfunc,...
%     covfunc,likfunc, x,y,xs);
% f15 = [mu15+2*sqrt(s215); flipdim(mu15-2*sqrt(s215),1)];
% f16 = [mu16+2*sqrt(s216); flipdim(mu16-2*sqrt(s216),1)];
% % ---------------------------------------------- %
% hyp17 = struct('mean',[], 'cov',[-10,0],'lik',-5);
% nlml17 = gp(hyp17, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% [mu17 s217] = gp(hyp17, @infGaussLik,meanfunc,...
%     covfunc,likfunc, x,y,xs);
% 
% hyp18 = minimize(hyp17, @gp, -100, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% nlml18 = gp(hyp18, @infGaussLik, meanfunc,...
%     covfunc, likfunc, x,y);
% [mu18 s218] = gp(hyp18, @infGaussLik,meanfunc,...
%     covfunc,likfunc, x,y,xs);
% f17 = [mu17+2*sqrt(s217); flipdim(mu17-2*sqrt(s217),1)];
% f18 = [mu18+2*sqrt(s218); flipdim(mu18-2*sqrt(s218),1)];
% 
% % ---------------------------------------------- %
% % Plot ----------------------------------------- %
% subplot(3,3,7);
% fill([xs; flipdim(xs,1)], f13, [5 5 5]/8)
% hold on 
% fill([xs; flipdim(xs,1)], f14, [7 7 7]/8)
% plot(xs, mu13); 
% plot(xs, mu14); 
% axis([-4 4 -3 4])
% plot(x, y, '+');
% title('[-0.25, 0],5');
% xlabel('input, x');
% ylabel('output, y');
% hold off
% % Plot ----------------------------------------- %
% subplot(3,3,8);
% fill([xs; flipdim(xs,1)], f15, [5 5 5]/8)
% hold on 
% fill([xs; flipdim(xs,1)], f16, [7 7 7]/8)
% plot(xs, mu15); 
% plot(xs, mu16); 
% axis([-4 4 -3 4])
% plot(x, y, '+');
% title('[-10, 0],-5');
% xlabel('input, x');
% ylabel('output, y');
% hold off
% % Plot ----------------------------------------- %
% subplot(3,3,9);
% fill([xs; flipdim(xs,1)], f17, [5 5 5]/8)
% hold on 
% fill([xs; flipdim(xs,1)], f18, [7 7 7]/8)
% plot(xs, mu17); 
% plot(xs, mu18); 
% axis([-4 4 -3 4])
% plot(x, y, '+');
% title('[10, 0],-5');
% xlabel('input, x');
% ylabel('output, y');
% hold off
% % Plot ----------------------------------------- %
% subplot(3,3,4);
% fill([xs; flipdim(xs,1)], f7, [5 5 5]/8)
% hold on 
% fill([xs; flipdim(xs,1)], f8, [7 7 7]/8)
% plot(xs, mu7); 
% plot(xs, mu8); 
% axis([-4 4 -3 4])
% plot(x, y, '+');
% title('[-0.25, -5],0');
% xlabel('input, x');
% ylabel('output, y');
% hold off
% % Plot ----------------------------------------- %
% subplot(3,3,5);
% fill([xs; flipdim(xs,1)], f9, [5 5 5]/8)
% hold on 
% fill([xs; flipdim(xs,1)], f10, [7 7 7]/8)
% plot(xs, mu9); 
% plot(xs, mu10); 
% axis([-4 4 -3 4])
% plot(x, y, '+');
% title('[-10, 5],0');
% xlabel('input, x');
% ylabel('output, y');
% hold off
% % Plot ----------------------------------------- %
% subplot(3,3,6);
% fill([xs; flipdim(xs,1)], f11, [5 5 5]/8)
% hold on 
% fill([xs; flipdim(xs,1)], f12, [7 7 7]/8)
% plot(xs, mu11); 
% plot(xs, mu12); 
% axis([-4 4 -3 4])
% plot(x, y, '+');
% title('[10, -5],0');
% xlabel('input, x');
% ylabel('output, y');
% hold off
% ---------------------------------------------- %
% Plot ----------------------------------------- %
% subplot(3,3,1);
% % fill([xs; flipdim(xs,1)], f0, [5 5 5]/8)
% % hold on 
% % fill([xs; flipdim(xs,1)], f, [7 7 7]/8)
% % plot(xs, mu0); 
% % plot(xs, mu); 
% % axis([-4 4 -3 4])
% % plot(x, y, '+');
% % % title('[-0.25, 0],0');
% % legend({'pre-train predictive error bar',...
% %         'post-train predictive error bar',...
% %         'pre-train predictive mean',...
% %         'post-train predictive mean',...
% %         'training set'},'location','bestoutside')
% % xlabel('input, x');
% % ylabel('output, y');
% % hold off
% % Plot ----------------------------------------- %
% subplot(3,3,2);
% fill([xs; flipdim(xs,1)], f3, [5 5 5]/8)
% hold on 
% fill([xs; flipdim(xs,1)], f4, [7 7 7]/8)
% plot(xs, mu3); 
% plot(xs, mu4); 
% axis([-4 4 -3 4])nln
% plot(x, y, '+');
% title('[-10, 0],0');
% xlabel('input, x');
% ylabel('output, y');
% hold off
% % Plot ----------------------------------------- %
% subplot(3,3,3);
% fill([xs; flipdim(xs,1)], f5, [5 5 5]/8)
% hold on 
% fill([xs; flipdim(xs,1)], f6, [7 7 7]/8)
% plot(xs, mu5); 
% plot(xs, mu6); 
% axis([-4 4 -3 4])
% plot(x, y, '+');
% title('[10, 0],0');
% xlabel('input, x');
% ylabel('output, y');
% hold off
% Table ---------------------------------------- %
% negative log marginal likelihood
Quantities = {'nlml', 'std.n', 'std.sf', 'length-scale'};
Pre_train = [nlml, exp(hyp.lik),exp(hyp.cov(2)),...
            exp(hyp.cov(1))];
Post_train = [nlml2,exp(hyp2.lik),exp(hyp2.cov(2)),...
                exp(hyp2.cov(1))];
disp(Post_train);

% ----------------------------------------------- %
% ------------------------------------------------- %
% d
% ------------------------------------------------- %
figure(4);
xd = linspace(-5,5,200)';
covfuncd = {@covProd, {@covPeriodic,@covSEiso}};
meanfuncd = @meanLinear;
likfuncd = @likGauss; sn = 0.1; 
hypd = struct('mean',0.5,'cov',[-0.5;0;0;2;0],...
    'lik',log(sn));
% without 1e-6*eye(200), K is not positive definite;
K = feval(covfuncd{:},hypd.cov,xd) + 1e-6*eye(200);
mud = feval(meanfuncd,hypd.mean,xd);
yd = chol(K)'*gpml_randn(0.15, 200, 1) + mud + exp(hypd.lik)*gpml_randn(0.2, 200, 1);
hypd2 = minimize(hypd, @gp, -100, @infGaussLik, meanfuncd,...
    covfuncd, likfuncd, xd,yd);
[mud s2d] = gp(hypd2, @infGaussLik,meanfuncd,...
    covfuncd,likfuncd, xd,yd,xs);
fd = [mud+2*sqrt(s2d); flipdim(mud-2*sqrt(s2d),1)];
fill([xs; flipdim(xs,1)], fd, [7 7 7]/8)
hold on 
plot(xs, mud);
plot(xd, yd, '+')
axis([-4 4 -6 6]);
legend({'post-train predictive error bar',...
        'post-train predictive mean',...
        'training set'},'location','bestoutside')
xlabel('input, x'); ylabel('output, y')
hold off

