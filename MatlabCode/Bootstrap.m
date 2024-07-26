clear;

reps = 1000;

% covids controls how bootstrap deals with 2020
% covids = 0: use pre-covid specification
% covids = 1: bootstrap from pre-covid, then add 2020 to dataset
% covids = 2: bootstrap from sample that includes covid
covids = 2;
corr_target = 0.175;
use_ridge = 1;

close all;
addpath("../ExternalCode/");

%use fixed random number seed for reproducibility
rng(12302018,'twister');

if covids > 0
    load("../Output/MainRun_Main.mat");
    Tmax = find(dts2>datetime(2019,12,31),1)-1;
    testSigma = opts.SigmaInit;
else
    load("../Output/MainRun_NoCOVID.mat");
    Tmax = T;
    testSigma = opts.SigmaInit;
end

sigs2 = sigs;
if opts.sigma_max <= 25
    sigs2 = [sigs2, opts.sigma_max+2:2:25];
    isigs2 = 1 ./ sigs2;
else
    sigs2 = sigs;
    isigs2 = 1 ./ sigs2;
end


cpi_full = table2array(Instruments(:, 'CPI')); 

RmZ = [Rm(:,2:T+2); Z(:,2:T+2); cpi_full(2:T+2)'];

Phi=([RmZ(:,1:T)',ones(T,1)] \ RmZ(:,2:T+1)');

epsRmZ =RmZ(:,2:T+1) - Phi'*[RmZ(:,1:T)',ones(T,1)]';

NRmZ = size(epsRmZ,1);
Nf = size(Rm,1);
Nr = size(R,1);

%analysis under null of zb rate = Rf

X = Rinput - iotaN*Rbinput; % excess returns of the portfolios
Xm = Rminput - iotaM*Rbinput; % excess returns of the factors
Xv = [Xm;ones(1,T)]';

BetaFull = Xv \ X';

epsR = X - BetaFull'*Xv';

%need to compute before changing gamma
wald = gamma' * inv(pvar(2:end-1,2:end-1)) * gamma;

if use_ridge == 1
    if covids > 0
        load("../Output/MainRun_Ridge.mat","ba","fitinfoa","zbrateReal","p_cons","gamma");
    else
        load("../Output/MainRun_RidgeNoCOVID.mat","ba","fitinfoa","zbrateReal","p_cons","gamma");
    end
    cbetas = [ba(:,fitinfoa.IndexMinMSE); fitinfoa.Intercept(fitinfoa.IndexMinMSE) ]
    p_cons = p_cons';
end

% epsC: residuals at point estimate
% epsC2: residuals if IID
epsC = cons_gr_ann' - cbetas'*[Zinput;ones(1,T)];

meanConsGr = mean(cons_gr_ann);
epsC2 = cons_gr_ann' - meanConsGr*ones(1,T);


%analysis under null of reduced correlation equation

%idea: gamma is as close a possible to point estimate with correlation
% equal to empirical real T-bill/consumption correlation
% reminder: in code, gamma is spread, not level (different from paper)
% also reminder: beta_inf is pre-normalization, and predicting -pi not pi
SigZ = Zinput * Zinput'/(T-1);
gamma_zbreal = gamma+ [1/ZSDiag(2,2); zeros(K-1,1)] + inv(ZSDiag(2:end,2:end))*beta_inf(1:end-1)*100;
gammaf = @(lam) inv((eye(K)-lam(1)*SigZ))*( gamma_zbreal + lam(2)*SigZ*cbetas(1:end-1));
%target = corr(RbinputReal',p_cons') * sqrt(gamma_zbreal'*SigZ*gamma_zbreal) * sqrt(cbetas(1:end-1)'*SigZ*cbetas(1:end-1));
target = corr_target * sqrt(gamma_zbreal'*SigZ*gamma_zbreal) * sqrt(cbetas(1:end-1)'*SigZ*cbetas(1:end-1));
consfunc = @(lam) [gammaf(lam)'*SigZ*gammaf(lam) - gamma_zbreal'*SigZ*gamma_zbreal, gammaf(lam)'*SigZ*cbetas(1:end-1) - target];

lamsol = fsolve(consfunc,[0,0]);

meanZbReal = mean(zbrateReal);
zbImpliedReal = [ones(1,T);Zinput]' * [meanZbReal; gammaf(lamsol)];
zbImplied = zbImpliedReal' - 100*(mean(exp_inf) - 1) - (Zinput'*inv(ZSDiag(2:end,2:end))*beta_inf(1:end-1)*100)';

%check
disp('targeted corr, std')
%corr(RbinputReal',p_cons')
corr_target
std(zbrateReal)
disp('calibration: ')
corr(zbImpliedReal,p_cons')
std(zbImpliedReal)

% run regressions under alternative hypothesis zero-beta rate
X2 = Rinput - iotaN*zbImplied; % excess returns of the portfolios
Xm2 = Rminput - iotaM*zbImplied; % excess returns of the factors
Xv2 = [Xm2;ones(1,T)]';

BetaFull2 = Xv2 \ X2';

epsR2 = X2 - BetaFull2'*Xv2';



eps = [epsRmZ; epsR; epsC; epsC2; epsR2];

walds = zeros(reps,1);
walds2 = zeros(reps,1);

Nsigs = length(sigs2);
svs2 = zeros(reps,Nsigs);

cFs = zeros(reps,1);
cFs2 = zeros(reps,1);

%some summary stats to check the bootstrap is working
meansRmZ = zeros(NRmZ,reps);
sdsRmZ = zeros(NRmZ,reps);
meansCons = zeros(2,reps);
sdsCons = zeros(2,reps);

meansR = zeros(Nr,reps);
sdsR = zeros(Nr,reps);

meansR2 = zeros(Nr,reps);
sdsR2 = zeros(Nr,reps);

ppool = gcp;
fevals(1:reps) = parallel.FevalFuture;

fhandle = @(ind) runSample (eps, T, Tmax, NRmZ, Phi, opts, iotaM, iotaN, BetaFull, ...
    meanConsGr, meanZbReal, lamsol, exp_inf, beta_inf, BetaFull2, ...
    testSigma, covids, sigs2, Nsigs, Nr, RmZ, Nf, cbetas, gammaf, ZSDiag, K);

for idx = 1:reps
    fevals(idx) = parfeval(ppool,fhandle,13,idx);
end

% hfig = waitbar(0,'Waiting...');
% 
% updateWaitbar = @(~) waitbar(mean({fevals.State} == "finished"),hfig);
% 
% updateWaitbarFutures = afterEach(fevals,updateWaitbar,0);
% afterAll(fevals,@(~) delete(hfig),0);

tic;
for idx = 1:reps
    [i,waldi, wald2, cF, cF2, sv2, meanRmZ, sdRmZ, meanCons, ...
     sdCons, meanR, sdR, meanR2, sdR2] = fetchNext(fevals);

    walds(i) = waldi;
    walds2(i) = wald2;
    cFs(i) = cF;
    cFs2(i) = cF2;
    svs2(i,:) = sv2;

    meansRmZ(:,i) = meanRmZ;
    sdsRmZ(:,i) = sdRmZ;
    meansCons(:,i) = meanCons;
    sdsCons(:,i) = sdCons;

    meansR(:,i) = meanR;
    meansR2(:,i) = meanR2;
    sdsR(:,i) = sdR;
    sdsR2(:,i) = sdR2;

    fprintf("Got result with index: %d.\n",i)
    toc;
end

% parfor i=1:reps
% 
%     [waldi, wald2, cF, cF2, sv2, meanRmZ, sdRmZ, meanCons, ...
%     sdCons, meanR, sdR, meanR2, sdR2] = runSample (eps, T, Tmax, NRmZ, Phi, opts, iotaM, iotaN, BetaFull, ...
%     meanConsGr, meanZbReal, lamsol, exp_inf, beta_inf, BetaFull2, ...
%     testSigma, covids, sigs2, Nsigs, Nr, RmZ, Nf, cbetas, gammaf, ZSDiag, K);
% 
%     walds(i) = waldi;
%     walds2(i) = wald2;
%     cFs(i) = cF;
%     cFs2(i) = cF2;
%     svs2(i,:) = sv2;
% 
%     meansRmZ(:,i) = meanRmZ;
%     sdsRmZ(:,i) = sdRmZ;
%     meansCons(:,i) = meanCons;
%     sdsCons(:,i) = sdCons;
% 
%     meansR(:,i) = meanR;
%     meansR2(:,i) = meanR2;
%     sdsR(:,i) = sdR;
%     sdsR2(:,i) = sdR2;
% 
% end

%print some checks
disp('mean and sd errors of Rm, Z, inflation')
mean(meansRmZ,2)-mean(RmZ,2)
mean(sdsRmZ,2) - std(RmZ,0,2)

disp('mean and sd errors of consumption growth')
mean(meansCons,2) - mean(cons_gr_ann)
mean(sdsCons,2) - std(cons_gr_ann)


save("../Output/BootstrapData-"+sprintf("%d",covids)+".mat");

pw = sum(walds>wald)/reps;

cfig=figure(1);
colororder({colors_list{2},colors_list{1}});
histogram(walds,'Normalization','pdf');
ylim([0,0.2]);
title('Bootstrap PDF of Wald Statistic under null of ZB=Rf+Cons.')
xline(wald, 'k-',{'Point Estimate p='+sprintf("%0.4f",pw)},'HandleVisibility','off');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapWald-"+sprintf("%d",covids)+".png");

mdl=fitlm(Zinput',cons_gr_ann);

cfig=figure(2);
colororder({colors_list{2},colors_list{1}});
histogram(cFs,'Normalization','pdf');
title('Bootstrap PDF of F Statistic under null of Cons. Unpredictable')
xline(mdl.ModelFitVsNullModel.Fstat, 'k-',{'Point Estimate p='+sprintf("%0.4f",mdl.ModelFitVsNullModel.Pvalue)},'HandleVisibility','off');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapF-"+sprintf("%d",covids)+".png");

ps2 = sum(svs2(:,sigs2==testSigma)>Svs(sigs2==testSigma))/reps;
pt2 = sum(svs2(:,sigs2==testSigma)>thresh)/reps;

cfig=figure(3);
colororder({colors_list{2},colors_list{1}});
histogram(svs2(:,sigs2==testSigma),'Normalization','pdf','BinWidth',thresh/4);
title('Bootstrap PDF of S-stat under alt. hypothesis of low correlation')
xline(Svs(sigs2==testSigma), 'k-',{'Point Estimate p='+sprintf("%0.4f",ps2)},'HandleVisibility','off');
xline(thresh, 'k-',{'Euler Rejection Threshold p='+sprintf("%0.4f",pt2)},'HandleVisibility','off');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapSstat2-"+sprintf("%d",covids)+".png");

cfig=figure(4);
colororder({colors_list{2},colors_list{1}});


xconf = [isigs2 isigs2(end:-1:1)] ;         
yconf = [log(prctile(svs2,95,1)) log(prctile(svs2(:,end:-1:1),5,1))];
p = fill(xconf,yconf,'c');

hold on;
title('Bootstrap of S-stat under alt. hypothesis of low correlation')
plot(isigs,log(Svs),'-.','LineWidth',2,'Color',colors_list{1});
plot(isigs2,log(mean(svs2,1)),'-','LineWidth',2,'Color',colors_list{2});
%plot(isigs,log(prctile(svs2,5,1)),'--','LineWidth',2,'Color',colors_list{2});
%plot(isigs,log(prctile(svs2,95,1)),'--','LineWidth',2,'Color',colors_list{2});
plot(isigs,log(threshs),'k:','LineWidth',1.5);
hold off;
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapSstat-"+sprintf("%d",covids)+".png");

pw2 = sum(walds2>wald)/reps;

cfig=figure(5);
colororder({colors_list{2},colors_list{1}});
histogram(walds2,'Normalization','pdf');
ylim([0,0.2]);
title('Bootstrap PDF of Wald Statistic under alt. hypothesis')
xline(wald, 'k-',{'Point Estimate'},'HandleVisibility','off');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapWald_alt-"+sprintf("%d",covids)+".png");

cfig=figure(6);
colororder({colors_list{2},colors_list{1}});
histogram(cFs2,'Normalization','pdf');
title('Bootstrap PDF of F Statistic under alt. hypothesis')
xline(mdl.ModelFitVsNullModel.Fstat, 'k-',{'Point Estimate'},'HandleVisibility','off');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapF_alt-"+sprintf("%d",covids)+".png");

pv_alt = sum(svs2>mean(thresh),1)/reps;
%pv_point = sum(svs2>Svs,1)/reps;

cfig=figure(7);
title('Bootstrap of S-stat rejecting prob. under alt. hypothesis of low correlation')
plot(isigs2,pv_alt,'-.','LineWidth',2,'Color',colors_list{1});
ylim([0,1]);
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapPower-"+sprintf("%d",covids)+".png");


function [wald, wald2, cF, cF2, sv2, meanRmZ, sdRmZ, meanCons, ...
    sdCons, meanR, sdR, meanR2, sdR2] = ...
    runSample (eps, T, Tmax, NRmZ, Phi, opts, iotaM, iotaN, BetaFull, ...
    meanConsGr, meanZbReal, lamsol, exp_inf, beta_inf, BetaFull2, ...
    testSigma, covids, sigs2, Nsigs, Nr, RmZ, Nf, cbetas, gammaf, ZSDiag, K)

 %draw new epsilons with replacement
    if covids == 1
        eps_sample = eps;
        eps_sample(:,1:Tmax) = datasample(eps(:,1:Tmax),Tmax,2);
    else
        eps_sample = datasample(eps,T,2);
    end

    %reconstruct time series
    %unpack residuals
    epsRmZ_sample = eps_sample(1:NRmZ,:);
    epsR_sample = eps_sample(NRmZ+1:NRmZ+Nr,:);
    epsC_sample = eps_sample(NRmZ+Nr+1,:);
    epsC2_sample = eps_sample(NRmZ+Nr+2,:);
    epsR2_sample = eps_sample(NRmZ+Nr+3:end,:);

    %run VAR forward
    RmZ_sample = zeros(size(RmZ));
    RmZ_sample(:,1) = RmZ(:,1);
    for j=1:T
        RmZ_sample(:,j+1) = epsRmZ_sample(:,j) + Phi'*[RmZ_sample(:,j)',1]';
    end

    %create sample
    Rminput_sample = RmZ_sample(1:Nf,2:end);
    Z_sample = RmZ_sample(Nf+1:end-1,1:end-1);
    cpi_sample = RmZ_sample(end,2:end);
    inflation_sample = log(1+cpi_sample/100);


    Rfex_sample = Z_sample(find(opts.Instruments=="RF"),1:end);
    Rbinput_sample = Rfex_sample;

    Xm_sample = Rminput_sample - iotaM*Rbinput_sample; % excess returns of the factors
    Xv_sample = [Xm_sample;ones(1,T)]';

    X_sample = epsR_sample + BetaFull'*Xv_sample';

    Rinput_sample = X_sample + iotaN*Rbinput_sample;

    %redo normalization 
    Zscales_sample = std(Z_sample');
    ZSDiag_sample = eye(length(Zscales_sample)+1);
    ZSDiag_sample(2:end,2:end) = diag(1./Zscales_sample);
    Zinput_sample = normalize(Z_sample')';
    
    cons_gr_ann_sample = epsC_sample + cbetas'*[Zinput_sample;ones(1,T)];
    cons_gr_ann_sample2 = epsC2_sample+meanConsGr*ones(1,T);

    cons_gr_ann_sample = cons_gr_ann_sample';
    cons_gr_ann_sample2 = cons_gr_ann_sample2';

    zbImpliedReal_sample = [ones(1,T);Zinput_sample]' * [meanZbReal; gammaf(lamsol)];
    zbImplied_sample = zbImpliedReal_sample' - 100*(mean(exp_inf) - 1) - (Zinput_sample'*inv(ZSDiag(2:end,2:end))*beta_inf(1:end-1)*100)';
    
    Xm2_sample = Rminput_sample - iotaM*zbImplied_sample; % excess returns of the factors
    Xv2_sample = [Xm2_sample;ones(1,T)]';

    X2_sample = epsR2_sample + BetaFull2'*Xv2_sample';

    Rinput2_sample = X2_sample + iotaN*zbImplied_sample;

    %construct F-stat to test consumption
    mdl_sample=fitlm(Zinput_sample',cons_gr_ann_sample2);
    cF=mdl_sample.ModelFitVsNullModel.Fstat;

    mdl_sample2=fitlm(Zinput_sample',cons_gr_ann_sample);
    cF2=mdl_sample2.ModelFitVsNullModel.Fstat;

    % reconstruct zero-beta rate under null of zbrate = rf + cons
    [~, Theta_sample, ~, ~, ~, ~, ~, ~,~,~,~,~,~,pvar_sample] = InstrumentGMMWrapperConc(Rinput_sample, Rminput_sample, Zinput_sample, Rbinput_sample, iotaN, iotaM, cons_gr_ann_sample/12, inflation_sample, Rfex_sample, 0,testSigma,'ZB',opts.har,opts.NLConsFactor,opts.SigmaType);

    gamma_sample = reshape(Theta_sample(2:1+K), K, 1);
    wald = gamma_sample' * inv(pvar_sample(2:end-1,2:end-1)) * gamma_sample;

    % reconstruct zero-beta rate S-test under null of low corr
    sv2 = zeros(1,Nsigs);
    for k = 1:Nsigs
    
        [~, Theta_sample_k, Sv_sample, ~, ~, ~, ~, ~,~,~,~,~,~,pvar_sample_k] = InstrumentGMMWrapperConc(Rinput2_sample, Rminput_sample, Zinput_sample, Rbinput_sample, iotaN, iotaM, cons_gr_ann_sample/12, inflation_sample, Rfex_sample, 0,sigs2(k),'ZB',opts.har,opts.NLConsFactor,opts.SigmaType);
        sv2(k) = Sv_sample;

        if sigs2(k) == testSigma
            gamma_sample2 = reshape(Theta_sample_k(2:1+K), K, 1);
            wald2 = gamma_sample2' * inv(pvar_sample_k(2:end-1,2:end-1)) * gamma_sample2;
        end
       
    end

    meanRmZ = mean(RmZ_sample,2);
    sdRmZ = std(RmZ_sample,0,2);
    meanCons = mean([cons_gr_ann_sample,cons_gr_ann_sample2]',2);
    sdCons = std([cons_gr_ann_sample,cons_gr_ann_sample2]',0,2);

    meanR = mean(Rinput_sample,2);
    meanR2 = mean(Rinput2_sample,2);
    sdR = std(Rinput_sample,0,2);
    sdR2 = std(Rinput2_sample,0,2);


end



