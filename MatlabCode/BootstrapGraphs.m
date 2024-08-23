
clear;
close all;
addpath("../ExternalCode/");
load("../Output/BootstrapData.mat");

pw = sum(walds(1:reps)>wald)/reps;

cfig=figure(1);
colororder({colors_list{2},colors_list{1}});
histogram(walds(1:reps),'Normalization','pdf');
ylim([0,0.2]);
xlim([0, 40]);
title('Bootstrap PDF of Wald Statistic under null of ZB=Rf+Cons.')
xline(wald, 'k-',{'Point Estimate p='+sprintf("%0.4f",pw)},'HandleVisibility','off');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapWald.png");


pwNC = sum(walds(reps+1:end)>waldNC)/reps;

cfig=figure(2);
colororder({colors_list{2},colors_list{1}});
hold on;
histogram(walds(reps+1:end),'Normalization','pdf');
ylim([0,0.2]);
xlim([0, 40]);
title('Bootstrap PDF of Wald Statistic under null of ZB=Rf+Cons.')
xline(wald, 'k-',{'Point Estimate p='+sprintf("%0.4f",pwNC)},'HandleVisibility','off');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapWaldNC.png");


%get the F effective
Zcov = cov(Zinput',0);
[coeffcov,~,coeff] = hac(Zinput',cons_gr_ann,Type="HC",Display="off");
cF = coeff(2:end)'*Zcov*coeff(2:end)/trace(coeffcov(2:end,2:end)*Zcov);

ZcovNC = cov(ZinputNC',0);
[coeffcovNC,~,coeffNC] = hac(ZinputNC',cons_gr_annNC,Type="HC",Display="off");
cFNC = coeffNC(2:end)'*ZcovNC*coeffNC(2:end)/trace(coeffcovNC(2:end,2:end)*ZcovNC);

pF = sum(cFs(1:reps)>cF)/reps;
pFNC = sum(cFs(reps+1:end)>cFNC)/reps;

%standard F stat
%mdl=fitlm(Zinput',cons_gr_ann);
%cF = mdl.ModelFitVsNullModel.Fstat;

pvNames = ["Wald","WaldNC","F-eff","F-effNC"];
pvals = [pw,pwNC,pF,pFNC];

writetable(table(pvNames',pvals'),"../Output/BootstrapPvals.csv");

cfig=figure(3);
colororder({colors_list{2},colors_list{1}});
histogram(cFs(1:reps),'Normalization','pdf');
title('Bootstrap PDF of F Statistic under null of Cons. Unpredictable')
xline(cF, 'k-',{'Point Estimate p='+sprintf("%0.4f",pF)},'HandleVisibility','off');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapF.png");

cfig=figure(4);
colororder({colors_list{2},colors_list{1}});
histogram(cFs(reps+1:end),'Normalization','pdf');
title('Bootstrap PDF of F Statistic under null of Cons. Unpredictable')
xline(cFNC, 'k-',{'Point Estimate p='+sprintf("%0.4f",pFNC)},'HandleVisibility','off');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapFNC.png");


ps2 = sum(svs2(1:reps,sigs2==testSigma)>Svs(sigs2==testSigma))/reps;
pt2 = sum(svs2(1:reps,sigs2==testSigma)>thresh)/reps;

cfig=figure(5);
colororder({colors_list{2},colors_list{1}});
histogram(svs2(1:reps,sigs2==testSigma),'Normalization','pdf','BinWidth',thresh/4);
title('Bootstrap PDF of S-stat under alt. hypothesis of low correlation')
xline(Svs(sigs2==testSigma), 'k-',{'Point Estimate p='+sprintf("%0.4f",ps2)},'HandleVisibility','off');
%xline(thresh, 'k-',{'Euler Rejection Threshold p='+sprintf("%0.4f",pt2)},'HandleVisibility','off');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapSstat2.png");

ps2NC_test = sum(svs2(reps+1:end,sigs2==testSigma)>SvsNC(sigs2==testSigma))/reps;
pt2NC_test = sum(svs2(reps+1:end,sigs2==testSigma)>thresh)/reps;

cfig=figure(6);
colororder({colors_list{2},colors_list{1}});
histogram(svs2(reps+1:end,sigs2==testSigma),'Normalization','pdf','BinWidth',thresh/4);
title('Bootstrap PDF of S-stat under alt. hypothesis of low correlation')
xline(SvsNC(sigs2==testSigma), 'k-',{'Point Estimate p='+sprintf("%0.4f",ps2NC_test)},'HandleVisibility','off');
%xline(thresh, 'k-',{'Euler Rejection Threshold p='+sprintf("%0.4f",pt2NC_test)},'HandleVisibility','off');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapSstat2NC.png");


ps2 = sum(svs2(1:reps,:)>Svs)/reps;
ps2NC = sum(svs2(reps+1:end,:)>SvsNC)/reps;

pv_alt99 = sum(svs2(1:reps,:)>mean(chi2inv(0.99,5)),1)/reps;
pv_altNC99 = sum(svs2(reps+1:end,:)>mean(chi2inv(0.99,5)),1)/reps;

pv_alt80 = sum(svs2(1:reps,:)>mean(chi2inv(0.80,5)),1)/reps;
pv_altNC80 = sum(svs2(reps+1:end,:)>mean(chi2inv(0.80,5)),1)/reps;

pv_alt = sum(svs2(1:reps,:)>mean(thresh),1)/reps;
pv_altNC = sum(svs2(reps+1:end,:)>mean(thresh),1)/reps;
%pv_point = sum(svs2>Svs,1)/reps;

cfig=figure(7);
colororder({colors_list{2},colors_list{1}});
title('Bootstrap of S-stat rejection prob. under alt. hypothesis of corr='+sprintf("%0.2f",corr_target));
hold on;
plot(isigs2,pv_alt80,'--x','LineWidth',2,'Color',colors_list{1});
plot(isigs2,pv_altNC80,'-x','LineWidth',2,'Color',colors_list{2});
plot(isigs2,pv_alt,'-.','LineWidth',2,'Color',colors_list{1});
plot(isigs2,pv_altNC,'-','LineWidth',2,'Color',colors_list{2});
plot(isigs2,pv_alt99,'.','LineWidth',2,'Color',colors_list{1});
plot(isigs2,pv_altNC99,'--','LineWidth',2,'Color',colors_list{2});

xscale log;
xticks([min(isigs2) 1/5 1/2 1 max(isigs2)]);
set(gca,'XMinorTick','Off');
hold off;
xlabel('IES $(1/\sigma)$, log scale','Interpreter','Latex');
ylabel('Prob. of S-stat > Threshold')
ylim([0,1.1]);
tightfig(cfig);
legend('80%, Incl. 2020','80%, Excl. 2020','95%, Incl. 2020','95%, Excl. 2020','99%, Incl. 2020','99%, Excl. 2020');
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapPower.png");


pvals = 1-chi2cdf(Svs,K);
pvalsNC = 1-chi2cdf(SvsNC,K);

cfig=figure(8);
hold on;
plot(isigs,pvals,'--','LineWidth',2,'Color',colors_list{1});
plot(isigs,pvalsNC,'-','LineWidth',2,'Color',colors_list{1});
plot(isigs,0.05*ones(size(pvals)),'k:','LineWidth',1.5);
xscale log;
xticks([min(isigs) 1/5 1/2 1 max(isigs)]);
ylim([0,1]);
set(gca,'XMinorTick','Off');
set(gca,'TickLabelInterpreter','latex')
%legend('Zero-Beta','T-Bill','Market','Threshold','Interpreter','Latex');
legend('Zero-Beta (Incl. 2020)','Zero-Beta (Excl. 2020)','Threshold','Interpreter','Latex');
ylabel('S-stat p-value','Interpreter','Latex');
xlabel('IES $(1/\sigma)$, log scale','Interpreter','Latex');
hold off;
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/StestPval_" + opts.Name + ".png");
