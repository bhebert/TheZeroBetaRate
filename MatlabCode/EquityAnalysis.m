close all;

%0: with CAPE
%1: CAPE, but gamma from Main
%2: Main
spec = 0;

if spec == 0 | spec == 1

    load("../Output/MainRun_AltCAPE.mat");
else
    load("../Output/MainRun_Main.mat");
end

if spec == 1

    gammaCAPE = gamma;
    load("../Output/MainRun_Main.mat",'gamma');

    gammaMAIN = gamma;
    index = find(contains(opts.Instruments,'CAPE'));

    gamma = zeros(size(gammaCAPE));
    gamma(1:index-1) = gammaMAIN(1:index-1);
    gamma(index+1:end) = gammaMAIN(index:end);
end

%print out some stats for the discussion section
mask=dts>datetime(2013,12,21) & dts<datetime(2019,1,1);
mean(zbrate(mask) - Rbinput(mask))*12

%one possible way to calibrate rho
%don't want to use DP because of divs vs buybacks issue
1/(1+mean(1./Instruments{:,"CAPE"}))

%a different way to calibrate rho
%implicitly assumes consumption claim return is zero-beta rate (on avg.)
exp(mean(cons_gr_ann/12-zbrateReal')/100)^12


%monthly value of rho
rho = (0.94)^(1/12);

% assumed value of 1/IES
sigma = 5;

%allows analysis at different horizons. not used in paper mainly for
%simplicity.
h=1;

% VAR analysis
% Run a VAR on the Zinput values
% recall the Zinput have already been centered here, so the inclusion of the
% constant does only a little
PhiZ=([Zinput(:,1:end-h)',ones(T-h,1)] \ Zinput(:,1+h:end)');

%gamma for nominal and real tbill rate
gamma_tbill = zeros(size(gamma));
gamma_tbill(1) = Zscales(1);
gamma_tbill_real = gamma_tbill + 100*(beta_inf(1:end-1).*Zscales');


%we need coefficients for the real zb rate, not the nominal or spread
%the beta_inf was created before the Z was normalized, hence the Zscales
%beta_inf predicts (1/(1+pi)) approx -pi, hence the + instead of -
gamma_real = gamma + gamma_tbill + 100*(beta_inf(1:end-1).*Zscales');

% Get the variable portion of the loading vector for the NPV calculation
%divide by hundred important, can't be percent returns
gvec = inv(eye(K)- (rho^h)*PhiZ(1:K,1:K))*gamma_real/100*h*(1/sigma-1);




gvec_tbill = inv(eye(K)- (rho^h)*PhiZ(1:K,1:K))*gamma_tbill_real/100*h*(1/sigma-1);

%compute the covariance matrix
%recall the Zinput have already been centered here
Varmat = Zinput(:,1:end)*Zinput(:,1:end)'/T;

%calculate the standard deviation
sdnpv = sqrt(gvec'*Varmat*gvec)

sdnpv_tbill = sqrt(gvec_tbill'*Varmat*gvec_tbill)

ts = Zinput'*gvec;

ts_tbill = Zinput'*gvec_tbill;

%redundant, just a double-check.
std(ts,1)

cape_series = table2array(Instruments(2:end-1,"CAPE"));

cfig=figure;
hold on;
plot(dts,ts,'-b','LineWidth',2);
plot(dts,ts_tbill,':g','LineWidth',2);
ylim([-1,1]);
ylabel("Log Valuation Ratio (Centered)")
plot(dts,log(cape_series)-mean(log(cape_series)),'-.r','LineWidth',2);
hold off;
legend({'VAR-implied Log PD Ratio (zero-beta)','VAR-implied Log PD Ratio (T-bill)','Log CAPE Ratio'},'Interpreter','Latex');



if spec == 0
    addpath("../ExternalCode/");
    tightfig(cfig);
    set(cfig,'PaperOrientation','landscape');
    print(cfig, '-dpng', "../Output/ValuationGraph.png");
end

