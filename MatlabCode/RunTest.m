clearvars -except opts optsDefault;
clc;
close all;

addpath("../ExternalCode/");

% color list
% visually distinct colors via colorbrewer2.org
colors_list = {'#e41a1c','#377eb8', '#4daf4a','#984ea3','#ff7f00'};

%% Load Data
Instruments = readtable("../Input/Instruments.csv");

Factors = readtable("../Input/Factors_Nominal.csv");

Portfolios = readtable("../Input/"+opts.AssetFile);


if strcmp(opts.ConsumptionVar,'JW')

    ConsData = readtable("../Input/Consumption/jw_m.csv");
    Cons_mask = table2array(ConsData(:,1)) >= datetime(1973,1,1) & table2array(ConsData(:,1)) <= datetime(2020,12,31);
    ConsGrData = ConsData(Cons_mask,[1,5]);

elseif strcmp(opts.ConsumptionVar,'Hall')

    ConsData = readtable("../Input/Consumption/hall_m.csv");
    Cons_mask = table2array(ConsData(:,1)) >= datetime(1973,1,1) & table2array(ConsData(:,1)) <= datetime(2020,12,31);
    ConsGrData = ConsData(Cons_mask,[1,5]); %4 is not per capita, 5 is per capita

else
    error('opts.ConsumptionVar must be one of JW or Hall');
end




%inflation
%this is month over month inflation
cpi = table2array(Instruments(:, 'CPI')); 

%for instruments, use one additional lag to avoid measurement error issue
%note: first value will be ignored later
Instruments(2:end,'CPI') = Instruments(1:end-1,'CPI');
Instruments(2:end,'CPI_rolling') = Instruments(1:end-1,'CPI_rolling');


R = table2array(Portfolios(:, 2:end))'; 
Rm = table2array(Factors(:, opts.Factors))'; 
Z = table2array(Instruments(:, opts.Instruments))'; 


Rf = table2array(Instruments(:, 'RF')); 
dates = Factors{:, 'Date'};



cons_gr_data = table2array(ConsGrData(:,2:end));

%Lag everything appropriately
%data starts in Jan 1973
%but need to lag; first return to be predicted is March 1973
%first instrument to be used is Feb 1973 (or Jan 1973 for 2x lagged
%variables)
Rinput = R(:, 3:end);
Rminput = Rm(:, 3:end);
T = size(Rminput,2);
Zinput = [Z(:, 2:end-1)];



%add lagged consumption
if opts.LaggedCons
    Zinput = [Zinput; cons_gr_data(1:end-2)'];
end


%this is supposed to be 1/(1+pi)
inf_factor = ones(size(cpi(1:end-2))) ./ (1+cpi(3:end)/100);

%log inflation
inflation = log(1+cpi(3:end)/100)';

%ex-post nominal return (= nominal yield) on t-bill over next month 
% i.e. Feb 1973 value is holding period return for March 1973
Rfex = Rf(2:end-1)';


%T-bill yield used as benchmark for zero-beta rate
Rbinput = Rf(2:end-1)';

%annualized real consumption growth
cons_gr_ann = 100*12*log(1+cons_gr_data(3:end)/100); 

%all assets are unit investment
iotaN = ones(size(Rinput,1), 1);

%all factors are zero-investment except market
iotaM = strcmp(opts.Factors,'Mkt')';

if opts.LinearConsumption
    Rminput = [Rminput; cons_gr_ann'];
    iotaM = [iotaM;0];
end

%dates for expected returns
dts = dates(2:end-1)';

%dates for realized returns
dts2 = dates(3:end)';


% code to truncate sample
if opts.NoCOVID 

    Tmax = find(dts2==datetime(2019,12,31));
  
    Rinput = Rinput(:,1:Tmax);  
    Zinput = Zinput(:,1:Tmax);
    Rminput = Rminput(:,1:Tmax);
    Rbinput = Rbinput(1:Tmax);
    dts = dts(1:Tmax);
    dts2 = dts(1:Tmax);
    cpi = cpi(1:Tmax+2);
    cons_gr_ann = cons_gr_ann(1:Tmax);
    inflation = inflation(1:Tmax);
    inf_factor = inf_factor(1:Tmax);
    Rfex = Rfex(1:Tmax);
    T=Tmax;
end



%predict the inflation factor
beta_inf = [Zinput',ones(size(cpi(1:end-2)))] \ inf_factor;
exp_inf =  [Zinput',ones(size(cpi(1:end-2)))] * beta_inf;


%expected real return on the Treasury bill
RbinputReal = 100*(Rbinput/100 + exp_inf' - 1);

%normalized inputs to unit variance
Zinput_pre = Zinput';
Zscales = std(Zinput');
ZSDiag = eye(length(Zscales)+1);
ZSDiag(2:end,2:end) = diag(1./Zscales);
Zinput = normalize(Zinput')';   % standardize the Z input variables.




%code to interact factors and Z
if opts.VaryingBeta
    Rm2 = permute(repmat(Rminput,[1,1,size(Zinput,1)]),[2,1,3]);
    Z2 = permute(repmat(Zinput,1,1,size(Rminput,1)),[2,3,1]);
    RZ2 = reshape(Rm2 .* Z2,T,[]);
    RminputZ = [Rminput;RZ2'];
    Rminput = RminputZ;
    iotaM = repmat(iotaM,1+size(Zinput,1),1);
end


%Values of Sigma to use in S-test
sigs = [0.25:0.25:0.75,1,4/3,5/3,2:1:opts.sigma_max];


% code for Ridge cross-validation to select penalty psi
if opts.RunRidge

    psis = linspace(0,199,50);

    Folds = 10; 
    L = length(psis);
    results = zeros(L, 2);

    parfor i = 1:L
        psii = psis(i);
        error = K_Fold(Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, cons_gr_ann/12, inflation, Rfex, psii,opts.SigmaInit,'ZB',opts.har,opts.NLConsFactor,Folds);
        results(i, :) = [psii,mean(error)];
    end

    results

    cfig = figure(4);
    plot(results(:,1), sqrt(results(:,2)),'LineWidth',2);
    xlabel('$\psi$', 'Interpreter','latex');
    ylabel('Loss');
    tightfig(cfig);

    set(cfig,'PaperOrientation','landscape');
    print(cfig, '-dpng', "Output/PsiGraph_"+opts.Name+".png");

    [v,ind] = min(results(:,2));
    psi = results(ind,1);
else
    psi = 0;
end

%Estimate zero-beta rate

[test, Theta, Sv, portRet3, weight, Sigma, Beta, alphas,mmoments,cmoments,amoments,thresh,mvar,pvar] = InstrumentGMMWrapperConc(Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, cons_gr_ann/12, inflation, Rfex, psi,opts.SigmaInit,'ZB',opts.har,opts.NLConsFactor);


[test; Sv; Theta]

K = size(Zinput,1);
Rf = Theta(1);
gamma = reshape(Theta(2:1+K), K, 1);
sigma = Theta(end);

rho = Theta(end-1);


[Rf;gamma]
sqrt(diag(pvar))

if ~opts.RunRidge
    [Rf;gamma] ./ sqrt(diag(pvar(1:end-1,1:end-1)))
else
    %if Ridge, no standard errors. This matrix will be ignored.
    pvar = eye(length([Rf;gamma])+1);
end

%save re-scaled coefficient estimates and standard errors
write(table(ZSDiag*[Rf;gamma], ZSDiag*pvar(1:end-1,1:end-1)*ZSDiag), "Output/ses_" + opts.Name + ".csv");


zbrate = gamma'*Zinput + Rf + Rbinput; % zero-beta rate

zbrateReal = zbrate + exp_inf' - 1;

%check that the mean and covariances make sense
mean(portRet3-zbrate)
std(portRet3-zbrate)
(portRet3-zbrate)*(Rminput-iotaM*zbrate)'/T

([Zinput;ones(size(portRet3))]' \ (portRet3-Rbinput)')' - [gamma;Rf]'

%save time series
write(table(dts2', RbinputReal', zbrateReal', portRet3' - inflation'*100, cons_gr_ann/12), "Output/zero_beta_rate" + opts.Name + ".csv")

% regression with both variables
fitlm([RbinputReal;zbrateReal]',cons_gr_ann/12)

% regression with only one variable
fitlm([zbrateReal]',cons_gr_ann/12)

fitlm([RbinputReal]',cons_gr_ann/12)

%if ridge, also use ridge for predicting consumption
if opts.RunRidge
    %note: this runs elastic net, alpha=1 is lasso, alpha->0 is ridge
    [ba,fitinfoa] = lasso(Zinput',cons_gr_ann,'CV',10,'Alpha',0.00001,'NumLambda',100000,'LambdaRatio',0.00000001);
    fitinfoa
    ba(:,fitinfoa.IndexMinMSE)
    p_cons = Zinput' * ba(:,fitinfoa.IndexMinMSE) + fitinfoa.Intercept(fitinfoa.IndexMinMSE);

else
    cbetas = [Zinput;ones(1,T)]' \ cons_gr_ann;
    p_cons = cbetas' * [Zinput;ones(1,T)];

end


fitlm(Zinput',portRet3 - inflation*100)
fitlm(Zinput',Rminput(1,:) - inflation*100)

%generate graph to compare ridge and non-ridge
if opts.RunRidge

    [~, Theta0] = InstrumentGMMWrapperConc(Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, cons_gr_ann/12, inflation, Rfex, 0,opts.SigmaInit,'ZB',opts.har,opts.NLConsFactor);
    Rf0 = Theta0(1);
    gamma0 = reshape(Theta0(2:1+K), K, 1);

    zbrate0 = gamma0'*Zinput + Rf0 + Rbinput; % zero-beta rate

    cfig=figure(3);
    colororder({colors_list{1},colors_list{4},colors_list{3},colors_list{2}});
    hold on
    plot(dts, zbrate*12, '-.', 'DisplayName', 'Zero-Beta Rate (Ridge)','LineWidth',2,'Color',colors_list{1});
    plot(dts, zbrate0*12, '-.', 'DisplayName', 'Zero-Beta Rate (Non-Ridge)','LineWidth',2,'Color',colors_list{4});
    plot(dts, Rbinput*12, '-.', 'DisplayName', 'T-Bill Yield.','LineWidth',2,'Color',colors_list{3});
    ylabel('Annualized Nominal Rate (mean +/- 4 s.d.)')
    legend('show')
    xlabel('Time')
    tightfig(cfig);
    set(cfig,'PaperOrientation','landscape');
    print(cfig, '-dpng', "Output/NominalRates_"+opts.Name+".png");
end

% plot the results.
cfig=figure(1);
colororder({colors_list{1},colors_list{2}});
hold on
yyaxis left;
plot(dts, zbrateReal*12, '-', 'DisplayName', 'Real Zero-Beta Rate','LineWidth',2,'Color',colors_list{1});
ylabel('Annualized Rate (mean +/- 4 s.d.)')

ylim(mean(zbrateReal'*12)+4*std(zbrateReal'*12)*[-1,1]);

yyaxis right;


ylim(mean(p_cons')+4*std(p_cons')*[-1,1]);
ylabel('Expect Cons. Growth (mean +/- 4 s.d.)');
plot(dts, p_cons, '--', 'DisplayName', 'Exp. Cons. Gr.','LineWidth',2,'Color',colors_list{2})
hold off
yline(mean(p_cons'), 'k-','HandleVisibility','off');
legend('show')
xlabel('Time')

tightfig(cfig);

set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "Output/ZBvCons_"+opts.Name+".png");


% plot the results.
cfig=figure(2);
colororder({colors_list{3},colors_list{2}});
hold on
yyaxis left;

plot(dts, RbinputReal*12, '-.', 'DisplayName', 'Exp. Real T-Bill Ret.','LineWidth',2,'Color',colors_list{3});
ylabel('Annualized Rate (mean +/- 4 s.d.)')
ylim(mean(RbinputReal'*12)+4*std(RbinputReal'*12)*[-1,1]);
yyaxis right;
ylim(mean(p_cons')+4*std(p_cons')*[-1,1]);
ylabel('Expect Cons. Growth (mean +/- 4 s.d.)');
plot(dts, p_cons, '--', 'DisplayName', 'Exp. Cons. Gr.','LineWidth',2,'Color',colors_list{2})
hold off
yline(mean(p_cons'), 'k-','HandleVisibility','off');
legend('show')
xlabel('Time')
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "Output/TBillvCons_"+opts.Name+".png");


%can't run GMM tests with ridge
if ~opts.RunRidge

Svs = zeros(size(sigs));
SvsRF = zeros(size(sigs));
SvsMkt = zeros(size(sigs));
tests = zeros(size(sigs));
threshs = zeros(size(sigs));


parfor i = 1:length(sigs)
    [ti,~,svi,~, ~, ~, ~, ~,~, ~,~,threshi] = InstrumentGMMWrapperConc(Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, cons_gr_ann/12, inflation,Rfex,0,sigs(i),'ZB',opts.har,opts.NLConsFactor);
    [~,~,sviRF] = InstrumentGMMWrapperConc(Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, cons_gr_ann/12, inflation,Rfex,0,sigs(i),'RF',opts.har,opts.NLConsFactor);
    [~,~,sviMkt] = InstrumentGMMWrapperConc(Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, cons_gr_ann/12, inflation,Rfex,0,sigs(i),'Mkt',opts.har,opts.NLConsFactor);
    tests(i) = ti;
    Svs(i) = svi;
    threshs(i)=threshi;
    SvsRF(i)=sviRF;
    SvsMkt(i) = sviMkt;
end

sigs(tests==1)

isigs = 1./sigs;

cfig=figure(3);
hold on;
plot(isigs,log(Svs),'-','LineWidth',2,'Color',colors_list{1});
plot(isigs,log(SvsRF),'-.','LineWidth',2,'Color',colors_list{3});
plot(isigs,log(SvsMkt),'--','LineWidth',2,'Color',colors_list{4});
plot(isigs,log(threshs),'k:','LineWidth',2);
legend('Zero-Beta','Tsy','Mkt','Threshold');
ylabel('log(S-stat)');
xlabel('ies (1/sigma)');
tightfig(cfig);
set(cfig,'PaperOrientation','landscape');
print(cfig, '-dpng', "Output/Stest_" + opts.Name + ".png");


end

%placebo figure, only in main analysis
if strcmp(opts.Name,'Main')


    bond_betas = [Zinput;ones(1,T)]' \ (Rminput(6,:)+Rbinput)';
    p_bond = bond_betas' * [Zinput;ones(1,T)] + exp_inf' - 1;
    

    
    cfig=figure(5);
    colororder({colors_list{4},colors_list{2}});
    hold on
    yyaxis left;

    plot(dts, p_bond*12, '-.', 'DisplayName', 'Exp. Real Tsy Bond Ret.','LineWidth',2,'Color',colors_list{4});
    ylabel('Annualized Rate (mean +/- 4 s.d.)')
    ylim(mean(p_bond'*12)+4*std(p_bond'*12)*[-1,1]);
    yyaxis right;
    ylim(mean(p_cons')+4*std(p_cons')*[-1,1]);
    ylabel('Expect Cons. Growth (mean +/- 4 s.d.)');
    plot(dts, p_cons, '--', 'DisplayName', 'Exp. Cons. Gr.','LineWidth',2,'Color',colors_list{2})
    hold off
    yline(mean(p_cons'), 'k-','HandleVisibility','off');
    legend('show')
    xlabel('Time')
    tightfig(cfig);
    set(cfig,'PaperOrientation','landscape');
    print(cfig, '-dpng', "Output/BondvCons_"+opts.Name+".png");
    
    
end
    
clear cfig ans;
save("Output/MainRun_"+ opts.Name +".mat");

