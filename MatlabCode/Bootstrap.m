
reps = 50;

close all;
addpath("../ExternalCode/");

%use fixed random number seed for reproducibility
rng(31416,'twister');

load("../Output/MainRun_Main.mat");

RmZ = [Rm(:,2:end); Z(:,2:end); cpi(2:end)'];

Phi=([RmZ(:,1:end-1)',ones(T,1)] \ RmZ(:,2:end)');

epsRmZ =RmZ(:,2:end) - Phi'*[RmZ(:,1:end-1)',ones(T,1)]';

NRmZ = size(epsRmZ,1);
Nf = size(Rm,1);
Nr = size(R,1);

%analysis under null of zb rate = Rf

X = Rinput - iotaN*Rbinput; % excess returns of the portfolios
Xm = Rminput - iotaM*Rbinput; % excess returns of the factors
Xv = [Xm;ones(1,T)]';

BetaFull = Xv \ X';

epsR = X - BetaFull'*Xv';

epsC = cons_gr_ann' - cbetas'*[Zinput;ones(1,T)];



%analysis under null of Euler equation

sdf_unscaled = exp(-opts.SigmaInit*cons_gr_ann/1200).*inf_factor;
delta = 1 / mean(sdf_unscaled) / (1+mean(zbrate)/100);
sdf = sdf_unscaled * delta;

func = @(beta,x) 1 ./ (x*beta);

beta_init = [1+mean(Rbinput)/100; ZSDiag(2,2)/100; zeros(K-1,1)];
sdf_model = fitnlm([ones(1,T);Zinput]',sdf,func,beta_init)

eps_sdf = (sdf - predict(sdf_model,[ones(1,T);Zinput]'))';

zbImplied = 100*([ones(1,T);Zinput]' * sdf_model.Coefficients.Estimate-1)';

X2 = Rinput - iotaN*zbImplied; % excess returns of the portfolios
Xm2 = Rminput - iotaM*zbImplied; % excess returns of the factors
Xv2 = [Xm2;ones(1,T)]';

BetaFull2 = Xv2 \ X2';

epsR2 = X2 - BetaFull2'*Xv2';



eps = [epsRmZ; epsR; epsC; eps_sdf; epsR2];

walds = zeros(reps,1);
svs = zeros(reps,1);

walds2 = zeros(reps,1);
svs2 = zeros(reps,1);

meansRmZ = zeros(NRmZ,reps);
sdsRmZ = zeros(NRmZ,reps);
meansCons = zeros(2,reps);
sdsCons = zeros(2,reps);
meanTrueZbs = zeros(2,reps);
sdsTrueZbs = zeros(2,reps);

meansR = zeros(Nr,reps);
sdsR = zeros(Nr,reps);

meansR2 = zeros(Nr,reps);
sdsR2 = zeros(Nr,reps);

parfor i=1:reps

    %draw new epsilons with replacement
    eps_sample = datasample(eps,T,2);

    %reconstruct time series
    epsRmZ_sample = eps_sample(1:NRmZ,:);
    epsR_sample = eps_sample(NRmZ+1:NRmZ+Nr,:);
    epsC_sample = eps_sample(NRmZ+Nr+1,:);
    epssdf_sample = eps_sample(NRmZ+Nr+2,:);
    epsR2_sample = eps_sample(NRmZ+Nr+3:end,:);


    RmZ_sample = zeros(size(RmZ));
    RmZ_sample(:,1) = RmZ(:,1);
    for j=1:T
        RmZ_sample(:,j+1) = epsRmZ_sample(:,j) + Phi'*[RmZ_sample(:,j)',1]';
    end

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

    Zscales_sample = std(Z_sample');
    ZSDiag_sample = eye(length(Zscales_sample)+1);
    ZSDiag_sample(2:end,2:end) = diag(1./Zscales_sample);
    Zinput_sample = normalize(Z_sample')';

    cons_gr_ann_sample = epsC_sample - cbetas'*[Zinput_sample;ones(1,T)];
    cons_gr_ann_sample = cons_gr_ann_sample';

    sdf_sample = epssdf_sample' + predict(sdf_model,[ones(1,T);Zinput_sample]');

    zbImplied_sample = 100*([ones(1,T);Zinput_sample]' * sdf_model.Coefficients.Estimate-1)';

    Xm2_sample = Rminput_sample - iotaM*zbImplied_sample; % excess returns of the factors
    Xv2_sample = [Xm2_sample;ones(1,T)]';

    X2_sample = epsR2_sample + BetaFull2'*Xv2_sample';

    Rinput2_sample = X2_sample + iotaN*zbImplied_sample;

    cons_gr_ann2_sample = -1200*log(sdf_sample / delta ./ exp(-inflation_sample'))/opts.SigmaInit;


    % reconstruct zero-beta rate under null of zbrate = rf + cons
    [~, Theta_sample, Sv_sample, ~, ~, ~, ~, ~,~,~,~,~,~,pvar_sample] = InstrumentGMMWrapperConc(Rinput_sample, Rminput_sample, Zinput_sample, Rbinput_sample, iotaN, iotaM, cons_gr_ann_sample/12, inflation_sample, Rfex_sample, 0,opts.SigmaInit,'ZB',opts.har,opts.NLConsFactor);

    % reconstruct zero-beta rate under null of Euler
    [~, Theta_sample2, Sv_sample2, ~, ~, ~, ~, ~,~,~,~,~,~,pvar_sample2] = InstrumentGMMWrapperConc(Rinput2_sample, Rminput_sample, Zinput_sample, Rbinput_sample, iotaN, iotaM, cons_gr_ann2_sample/12, inflation_sample, Rfex_sample, 0,opts.SigmaInit,'ZB',opts.har,opts.NLConsFactor);

    Rf_sample = Theta_sample(1);
    gamma_sample = reshape(Theta_sample(2:1+K), K, 1);

    Rf_sample2 = Theta_sample2(1);
    gamma_sample2 = reshape(Theta_sample2(2:1+K), K, 1);

    walds(i) = gamma_sample' * inv(pvar_sample(2:end-1,2:end-1)) * gamma_sample;
    svs(i) = Sv_sample;

    walds2(i) = gamma_sample2' * inv(pvar_sample2(2:end-1,2:end-1)) * gamma_sample2;
    svs2(i) = Sv_sample2;

    meansRmZ(:,i) = mean(RmZ_sample,2);
    sdsRmZ(:,i) = std(RmZ_sample,0,2);
    meansCons(:,i) = mean([cons_gr_ann_sample,cons_gr_ann2_sample]',2);
    sdsCons(:,i) = std([cons_gr_ann_sample,cons_gr_ann2_sample]',0,2);
    meanTrueZbs(:,i) = mean([Rbinput_sample;zbImplied_sample],2);
    sdsTrueZbs(:,i) = std([Rbinput_sample;zbImplied_sample],0,2);

    meansR(:,i) = mean(Rinput_sample,2);
    meansR2(:,i) = mean(Rinput2_sample,2);
    sdsR(:,i) = std(Rinput_sample,0,2);
    sdsR2(:,i) = std(Rinput2_sample,0,2);
end

wald = gamma' * inv(pvar(2:end-1,2:end-1)) * gamma;

pw = sum(walds>wald)/reps;

cfig=figure(1);
colororder({colors_list{2},colors_list{1}});
histogram(walds,'Normalization','pdf');
ylim([0,0.2]);
title('Bootstrap PDF of Wald Statistic under null of ZB=Rf+Cons.')
xline(wald, 'k-',{'Point Estimate p='+sprintf("%0.4f",pw)},'HandleVisibility','off');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapWald.png");

ps = sum(svs<=Sv)/reps;
pt = sum(svs<=thresh)/reps;

cfig=figure(2);
colororder({colors_list{2},colors_list{1}});
histogram(svs,'Normalization','pdf','BinWidth',thresh/4);
ylim([0,0.2]);
xlim([0,50]);
title('Bootstrap PDF of S-stat under null of ZB=Rf+Cons.')
xline(Sv, 'k-',{'Point Estimate p='+sprintf("%0.4f",ps)},'HandleVisibility','off');
xline(thresh, 'k-',{'Euler Rejection Threshold p='+sprintf("%0.4f",pt)},'HandleVisibility','off');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapSstat.png");

pw2 = sum(walds2>wald)/reps;

cfig=figure(3);
colororder({colors_list{2},colors_list{1}});
histogram(walds2,'Normalization','pdf');
ylim([0,0.2]);
title('Bootstrap PDF of Wald Statistic under null of Euler.')
xline(wald, 'k-',{'Point Estimate p='+sprintf("%0.4f",pw2)},'HandleVisibility','off');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapWald2.png");

ps2 = sum(svs2<=Sv)/reps;
pt2 = sum(svs2<=thresh)/reps;

cfig=figure(4);
colororder({colors_list{2},colors_list{1}});
histogram(svs2,'Normalization','pdf','BinWidth',thresh/4);
ylim([0,0.2]);
xlim([0,50]);
title('Bootstrap PDF of S-stat under null of Euler.')
xline(Sv, 'k-',{'Point Estimate p='+sprintf("%0.4f",ps2)},'HandleVisibility','off');
xline(thresh, 'k-',{'Euler Rejection Threshold p='+sprintf("%0.4f",pt2)},'HandleVisibility','off');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapSstat2.png");

