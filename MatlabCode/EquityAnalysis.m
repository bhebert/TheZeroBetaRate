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


%monthly value of rho
rho = (0.96)^(1/12);

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


% Get the variable portion of the loading vector for the NPV calculation
%divide by hundred important, can't be percent returns
gvec = inv(eye(K)- (rho^h)*PhiZ(1:K,1:K))*gamma/100*h*(1/sigma-1);

%compute the covariance matrix
%recall the Zinput have already been centered here
Varmat = Zinput(:,1:end)*Zinput(:,1:end)'/T;

%calculate the standard deviation
sdnpv = sqrt(gvec'*Varmat*gvec)

ts = Zinput'*gvec;

%redundant, just a double-check.
std(ts,1)

cape_series = table2array(Instruments(2:end-1,"CAPE"));

cfig=figure;
hold on;
plot(dts,ts,'-b','LineWidth',2);
ylim([-1,1]);
ylabel("Log Valuation Ratio (Centered)")
plot(dts,log(cape_series)-mean(log(cape_series)),'-.r','LineWidth',2);
hold off;
legend({'VAR-implied Log PD Ratio of Consumption Claim','Log CAPE Ratio'},'Interpreter','Latex');



if spec == 0
    addpath("../ExternalCode/");
    tightfig(cfig);
    set(cfig,'PaperOrientation','landscape');
    print(cfig, '-dpng', "../Output/ValuationGraph.png");
end

