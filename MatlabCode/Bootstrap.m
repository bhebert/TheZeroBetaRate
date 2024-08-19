clear;

%for testing
reps = 1000;

use_ridge = 0;

close all;
addpath("../ExternalCode/");

%use fixed random number seed for reproducibility
rng(12302018,'twister');

load("../Output/MainRun_Main.mat");
testSigma = 8; %opts.SigmaInit;




sigs2 = sigs;
isigs2 = 1 ./ sigs2;
% if opts.sigma_max <= 25
%     sigs2 = [sigs2, opts.sigma_max+2:2:25];
%     isigs2 = 1 ./ sigs2;
% else
%     sigs2 = sigs;
%     isigs2 = 1 ./ sigs2;
% end

%the code shifts the Instruments array for CPI, so don't use this
%cpi_full = table2array(Instruments(:, 'CPI')); 
cpi_full = cpi;


%time series of VAR state variables (factors, instruments, inflation)
RmZ = [Rm(:,2:T+2); Z(:,2:T+2); cpi_full(2:T+2)'];

%run regression excluding covid
Phi=([RmZ(:,1:Tmax)',ones(Tmax,1)] \ RmZ(:,2:Tmax+1)');

%construct residuals for full and NC samples
epsRmZ =RmZ(:,2:T+1) - Phi'*[RmZ(:,1:T)',ones(T,1)]';
epsRmZNC =RmZ(:,2:Tmax+1) - Phi'*[RmZ(:,1:Tmax)',ones(Tmax,1)]';

NRmZ = size(epsRmZ,1);
Nf = size(Rm,1);
Nr = size(R,1);

%analysis under null of zb rate = Rf

X = Rinput - iotaN*Rbinput; % excess returns of the portfolios
Xm = Rminput - iotaM*Rbinput; % excess returns of the factors
Xv = [Xm;ones(1,T)]';

%compute betas from no-covid sample
BetaFull = Xv(1:Tmax,:) \ X(:,1:Tmax)';

%construct return residuals w/ and w/o covid
epsR = X - BetaFull'*Xv';
epsRNC = epsR(:,1:Tmax);

%compute wald statistics
%need to compute before changing gamma (if ridge gamma is used)
wald = gamma' * inv(pvar(2:end-1,2:end-1)) * gamma;
waldNC = gammaNC' * inv(pvarNC(2:end-1,2:end-1)) * gammaNC;

if use_ridge == 1

    load("../Output/MainRun_Ridge.mat","baNC","fitinfoaNC","zbrateRealNC","p_consNC","gammaNC","RbinputRealNC");
    
    %overwrite consumption betas
    cbetasNC = [baNC(:,fitinfoaNC.IndexMinMSE); fitinfoaNC.Intercept(fitinfoaNC.IndexMinMSE) ]
    p_consNC = p_consNC';
end

%choose correlation target
corr_target = corr(RbinputRealNC',p_consNC');

% epsC: residuals at point estimate
% epsC2: residuals if IID
epsC = cons_gr_ann' - cbetasNC'*[Zinput;ones(1,T)];
epsCNC = epsC(1:Tmax);

meanConsGr = mean(cons_gr_annNC);
epsC2 = cons_gr_ann' - meanConsGr*ones(1,T);
epsC2NC = epsC2(1:Tmax);

%analysis under null of reduced correlation equation

%idea: gamma is as close a possible to point estimate with correlation
% equal to empirical real T-bill/consumption correlation
% reminder: in code, gamma is spread, not level (different from paper)
% also reminder: beta_inf is pre-normalization, and predicting (approx.) -pi not pi
SigZ = ZinputNC * ZinputNC'/(Tmax-1);

%gamma for real zb rate
gamma_zbreal = gammaNC+ [1/ZSDiagNC(2,2); zeros(K-1,1)] + inv(ZSDiagNC(2:end,2:end))*beta_infNC(1:end-1)*100;

%first-order condition of the minimization problem, gamma (zb real) as a function of
%multipliers lam(1) and lam(2)
%updated FOC under new norm
gammaf = @(lam) 1/(1-lam(1))*( gamma_zbreal + lam(2)*cbetasNC(1:end-1));

%covariance target
target = corr_target * sqrt(gamma_zbreal'*SigZ*gamma_zbreal) * sqrt(cbetasNC(1:end-1)'*SigZ*cbetasNC(1:end-1));

%constraints in the minimization problem: match variance of real zb rate
%and targeted covariance with consumption
consfunc = @(lam) [gammaf(lam)'*SigZ*gammaf(lam) - gamma_zbreal'*SigZ*gamma_zbreal, gammaf(lam)'*SigZ*cbetasNC(1:end-1) - target];

%solve for multipliers
[lamsol,~,flag] = fsolve(consfunc,[0,0]);

assert(flag > 0,"Failed to solve for adjusted gamma parameters");

%reconstruct zb rate from rotated parameters
%note that the parameters are estimated on the NC sample
%but the reconstructed series has the full sample length
%because of this, we need to create ZinputFS, which is the NC Zinput + out
%of sample
ZinputFS = ZSDiagNC(2:end,2:end)*(Z(:, 2:end-1) - mean(Z(:,2:1+Tmax),2));
meanZbReal = mean(zbrateRealNC);
zbImpliedReal = [ones(1,T);ZinputFS]' * [meanZbReal; gammaf(lamsol)];
zbImplied = zbImpliedReal' - 100*(mean(exp_infNC) - 1) - (ZinputFS'*inv(ZSDiagNC(2:end,2:end))*beta_infNC(1:end-1)*100)';

%check
disp('targeted corr, std')
%corr(RbinputReal',p_cons')
corr_target
std(zbrateRealNC)
disp('calibration: ')
corr(zbImpliedReal(1:Tmax),p_consNC')
std(zbImpliedReal(1:Tmax))

% run regressions under alternative hypothesis zero-beta rate
X2 = Rinput - iotaN*zbImplied; % excess returns of the portfolios
Xm2 = Rminput - iotaM*zbImplied; % excess returns of the factors
Xv2 = [Xm2;ones(1,T)]';

BetaFull2 = Xv2(1:Tmax,:) \ X2(:,1:Tmax)';

epsR2 = X2 - BetaFull2'*Xv2';
epsR2NC = epsR2(:,1:Tmax);



eps = [epsRmZ; epsR; epsC; epsC2; epsR2];
epsNC = [epsRmZNC; epsRNC; epsCNC; epsC2NC; epsR2NC];

walds = zeros(2*reps,1);
walds2 = zeros(2*reps,1);



Nsigs = length(sigs2);
svs2 = zeros(2*reps,Nsigs);


cFs = zeros(2*reps,1);
cFs2 = zeros(2*reps,1);

%some summary stats to check the bootstrap is working
meansRmZ = zeros(NRmZ,2*reps);
sdsRmZ = zeros(NRmZ,2*reps);
meansCons = zeros(2,2*reps);
sdsCons = zeros(2,2*reps);

meansR = zeros(Nr,2*reps);
sdsR = zeros(Nr,2*reps);

meansR2 = zeros(Nr,2*reps);
sdsR2 = zeros(Nr,2*reps);

ppool = gcp;
fevals(1:2*reps) = parallel.FevalFuture;

fhandle = @(ind) runSample (eps, T, NRmZ, Phi, opts, iotaM, iotaN, BetaFull, ...
    meanConsGr, meanZbReal, lamsol, exp_infNC, beta_infNC, BetaFull2, ...
    testSigma, sigs2, Nsigs, Nr, RmZ, Nf, cbetasNC, gammaf, ZSDiagNC, K);

fhandleNC = @(ind) runSample (epsNC, Tmax, NRmZ, Phi, opts, iotaM, iotaN, BetaFull, ...
    meanConsGr, meanZbReal, lamsol, exp_infNC, beta_infNC, BetaFull2, ...
    testSigma, sigs2, Nsigs, Nr, RmZ(:,1:Tmax+1), Nf, cbetasNC, gammaf, ZSDiagNC, K);

%test a function call
%fhandleNC(0)

fprintf("Done with test")

for idx = 1:reps
    fevals(idx) = parfeval(ppool,fhandle,13,idx);
    fevals(idx+reps) = parfeval(ppool,fhandleNC,13,idx);
end

tic;
for idx = 1:2*reps
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
mean(meansRmZ(:,reps+1:end),2)-mean(RmZ(:,2:Tmax+1),2)
mean(sdsRmZ(:,reps+1:end),2) - std(RmZ(:,2:Tmax+1),0,2)

disp('mean and sd errors of consumption growth')
mean(meansCons(:,reps+1:end),2) - mean(cons_gr_annNC)
mean(sdsCons(:,reps+1:end),2) - std(cons_gr_annNC)

clear fevals;
save("../Output/BootstrapData.mat");




function [wald, wald2, cF, cF2, sv2, meanRmZ, sdRmZ, meanCons, ...
    sdCons, meanR, sdR, meanR2, sdR2] = ...
    runSample (eps, T, NRmZ, Phi, opts, iotaM, iotaN, BetaFull, ...
    meanConsGr, meanZbReal, lamsol, exp_infNC, beta_infNC, BetaFull2, ...
    testSigma, sigs2, Nsigs, Nr, RmZ, Nf, cbetasNC, gammaf, ZSDiagNC, K)

     %draw new epsilons with replacement
     eps_sample = datasample(eps,T,2);
    
    
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
    
    cons_gr_ann_sample = epsC_sample + cbetasNC'*[Zinput_sample;ones(1,T)];
    cons_gr_ann_sample2 = epsC2_sample+meanConsGr*ones(1,T);

    cons_gr_ann_sample = cons_gr_ann_sample';
    cons_gr_ann_sample2 = cons_gr_ann_sample2';

    zbImpliedReal_sample = [ones(1,T);Zinput_sample]' * [meanZbReal; gammaf(lamsol)];
    zbImplied_sample = zbImpliedReal_sample' - 100*(mean(exp_infNC) - 1) - (Zinput_sample'*inv(ZSDiagNC(2:end,2:end))*beta_infNC(1:end-1)*100)';
    
    Xm2_sample = Rminput_sample - iotaM*zbImplied_sample; % excess returns of the factors
    Xv2_sample = [Xm2_sample;ones(1,T)]';

    X2_sample = epsR2_sample + BetaFull2'*Xv2_sample';

    Rinput2_sample = X2_sample + iotaN*zbImplied_sample;

    %construct F-stat to test consumption
    %remember that sample2 for cons_gr is the IID assumption
    %while sample is the point estimate. this is why the code below
    % looks like it "switches" 1 and 2

    %non-robust F
    %mdl_sample=fitlm(Zinput_sample',cons_gr_ann_sample2);
    %cF=mdl_sample.ModelFitVsNullModel.Fstat;

    %mdl_sample2=fitlm(Zinput_sample',cons_gr_ann_sample);
    %cF2=mdl_sample2.ModelFitVsNullModel.Fstat;

    %effective F with robust std. errors
    Zcov = cov(Zinput_sample',0);
    [coeffcov,~,coeff] = hac(Zinput_sample',cons_gr_ann_sample2,Type="HC",Display="off");
    cF = coeff(2:end)'*Zcov*coeff(2:end)/trace(coeffcov(2:end,2:end)*Zcov);

    [coeffcov2,~,coeff2] = hac(Zinput_sample',cons_gr_ann_sample,Type="HC",Display="off");
    cF2 = coeff2(2:end)'*Zcov*coeff2(2:end)/trace(coeffcov2(2:end,2:end)*Zcov);

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



