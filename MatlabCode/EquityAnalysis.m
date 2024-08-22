close all;
clear;

%0: with CAPE
%1: CAPE, but gamma from Main
%2: Main
%3: with CAPE and Ridge
spec = 0;

if spec == 0 | spec == 1

    load("../Output/MainRun_AltCAPE.mat");
elseif spec == 3
    load("../Output/MainRun_AltCAPERidge.mat");
else
    load("../Output/MainRun_Main.mat");
end

if spec == 1

    gammaCAPE = gammaNC;
    load("../Output/MainRun_Main.mat",'gammaNC');

    gammaMAIN = gammaNC;
    index = find(contains(opts.Instruments,'CAPE'));

    gammaNC = zeros(size(gammaCAPE));
    gammaNC(1:index-1) = gammaMAIN(1:index-1);
    gammaNC(index+1:end) = gammaMAIN(index:end);
end

%print out some stats for the discussion section
mask=dtsNC>datetime(2013,12,21) & dtsNC<datetime(2019,1,1);
mean(zbrateNC(mask) - RbinputNC(mask))*12

%one possible way to calibrate rho
%don't want to use DP because of divs vs buybacks issue
1/(1+mean(1./Instruments{:,"CAPE"}))

%a different way to calibrate rho
%implicitly assumes consumption claim return is zero-beta rate (on avg.)
exp(mean(cons_gr_annNC/12-zbrateRealNC')/100)^12


%monthly value of rho
rho = (0.94)^(1/12);

% assumed value of 1/IES
sigma = 7.5;

%allows analysis at different horizons. not used in paper mainly for
%simplicity.
h=1;

% VAR analysis
% Run a VAR on the Zinput values
% recall the Zinput have already been centered here, so the inclusion of the
% constant does only a little
PhiZ=([ZinputNC(:,1:end-h)',ones(Tmax-h,1)] \ ZinputNC(:,1+h:end)');

%gamma for nominal and real tbill rate
gamma_tbill = zeros(size(gamma));
gamma_tbill(1) = ZscalesNC(1);
gamma_tbill_real = gamma_tbill + 100*(beta_infNC(1:end-1).*ZscalesNC');


%we need coefficients for the real zb rate, not the nominal or spread
%the beta_inf was created before the Z was normalized, hence the Zscales
%beta_inf predicts (1/(1+pi)) approx -pi, hence the + instead of -
gamma_real = gammaNC + gamma_tbill + 100*(beta_infNC(1:end-1).*ZscalesNC');

% Get the variable portion of the loading vector for the NPV calculation
%divide by hundred important, can't be percent returns
gvec = inv(eye(K)- (rho^h)*PhiZ(1:K,1:K))*gamma_real/100*h*(1/sigma-1);




gvec_tbill = inv(eye(K)- (rho^h)*PhiZ(1:K,1:K))*gamma_tbill_real/100*h*(1/sigma-1);

%compute the covariance matrix
%recall the Zinput have already been centered here
Varmat = ZinputNC(:,1:end)*ZinputNC(:,1:end)'/Tmax;

%calculate the standard deviation
sdnpv = sqrt(gvec'*Varmat*gvec)

sdnpv_tbill = sqrt(gvec_tbill'*Varmat*gvec_tbill)

ts = ZinputNC'*gvec;

ts_tbill = ZinputNC'*gvec_tbill;

%redundant, just a double-check.
std(ts,1)

cape_series = table2array(Instruments(2:end-1,"CAPE"));
cape_series = cape_series(1:Tmax);

save("../Output/EquityAnalysis_"+ opts.Name +".mat");

cfig=figure;
hold on;
plot(dtsNC,ts,'-b','LineWidth',2);
plot(dtsNC,ts_tbill,':g','LineWidth',2);
ylim([-1,1]);
ylabel("Log Valuation Ratio (Centered)")
plot(dtsNC,log(cape_series)-mean(log(cape_series)),'-.r','LineWidth',2);
hold off;
legend({'VAR-implied Log PD Ratio (zero-beta)','VAR-implied Log PD Ratio (T-bill)','Log CAPE Ratio'},'Interpreter','Latex');



if spec == 0
    addpath("../ExternalCode/");
    tightfig(cfig);
    set(cfig,'PaperOrientation','landscape');
    print(cfig, '-dpng', "../Output/ValuationGraph.png");
end

