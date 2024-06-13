
reps = 10;

close all;
addpath("../ExternalCode/");

load("../Output/MainRun_Main.mat");

RmZ = [Rm(:,2:end); Z(:,2:end)];

Phi=([RmZ(:,1:end-1)',ones(T,1)] \ RmZ(:,2:end)');

epsRmZ =RmZ(:,2:end) - Phi'*[RmZ(:,1:end-1)',ones(T,1)]';

NRmZ = size(epsRmZ,1);
Nf = size(Rm,1);

%analysis under null of zb rate = Rf

X = Rinput - iotaN*Rbinput; % excess returns of the portfolios
Xm = Rminput - iotaM*Rbinput; % excess returns of the factors
Xv = [Xm;ones(1,T)]';

BetaFull = Xv \ X';

epsR = X - BetaFull'*Xv';

epsC = cons_gr_ann' - cbetas'*[Zinput;ones(1,T)];

eps = [epsRmZ; epsR; epsC];

walds = zeros(reps,1);
svs = zeros(reps,1);

parfor i=1:reps

    %draw new epsilons with replacement
    eps_sample = datasample(eps,T,2);

    %reconstruct time series
    epsRmZ_sample = eps_sample(1:NRmZ,:);
    epsR_sample = eps_sample(NRmZ+1:end-1,:);
    epsC_sample = eps_sample(end,:);

    RmZ_sample = zeros(size(RmZ));
    RmZ_sample(:,1) = RmZ(:,1);
    for j=1:T
        RmZ_sample(:,j+1) = epsRmZ_sample(:,j) + Phi'*[RmZ_sample(:,j)',1]';
    end

    Rminput_sample = RmZ_sample(1:Nf,2:end);
    Z_sample = RmZ_sample(Nf+1:end,1:end-1);

    Rfex_sample = Z_sample(find(opts.Instruments=="RF"),1:end);
    Rbinput_sample = Rfex;

    %inflation only is used in specificationa with non-linear cons factor
    %so setting it to zero will not affect the results
    inflation_sample = zeros(size(inflation));

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

    % reconstruct zero-beta rate
    [~, Theta_sample, Sv_sample, ~, ~, ~, ~, ~,~,~,~,~,~,pvar_sample] = InstrumentGMMWrapperConc(Rinput_sample, Rminput_sample, Zinput_sample, Rbinput_sample, iotaN, iotaM, cons_gr_ann_sample/12, inflation_sample, Rfex_sample, 0,opts.SigmaInit,'ZB',opts.har,opts.NLConsFactor);

    Rf_sample = Theta_sample(1);
    gamma_sample = reshape(Theta_sample(2:1+K), K, 1);

    walds(i) = gamma_sample' * inv(pvar_sample(2:end-1,2:end-1)) * gamma_sample;
    svs(i) = Sv_sample;

end

wald = gamma' * inv(pvar(2:end-1,2:end-1)) * gamma;

cfig=figure(1);
colororder({colors_list{2},colors_list{1}});
histogram(walds,'Normalization','pdf');
ylim([0,0.2]);
title('Bootstrap PDF of Wald Statistic under null of ZB=Rf+Cons.')
xline(wald, 'k-',{'Point Estimate'},'HandleVisibility','off');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapWald.png");


cfig=figure(2);
colororder({colors_list{2},colors_list{1}});
histogram(svs,'Normalization','pdf','BinWidth',thresh/2);
ylim([0,0.2]);
xlim([0,50]);
title('Bootstrap PDF of S-stat under null of ZB=Rf+Cons.')
xline(Sv, 'k-',{'Point Estimate'},'HandleVisibility','off');
xline(thresh, 'k-',{'Euler Rejection Threshold'},'HandleVisibility','off');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "../Output/BootstrapSstat.png");

