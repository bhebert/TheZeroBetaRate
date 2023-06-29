%control code for generating results





%Name: name for file extensions, string

%AssetFile: name for Input/xxx to read asset returns
% 27_plus_Industry_Portfolios_Nominal.csv
% FF5_plus_Industry_Portfolios_Nominal.csv

%Factors: name of factor set to include. Options:
% Mkt, SMB, HML, RMW, CMA, term_spread, DEF

%NLConsFactor: true/false, included non-linear consumption factor

%VaryingBeta: true/false, add Instruments x Factors as factors

%LinearConsumption: true or false, include linear consumption as a factor
% if true, probably want SigmaInit = 0

%SigmaInit: Value of Sigma for non-linear consumption factor
% >0: run with this value
% 0: inflation only

%Instruments:
%Options are
% RF       CPI       UMP        EBP         CAPE      TSP     CPI_rolling    shadow_spread     BAAS


%LaggedCons: true/false, add lagged consumption growth as instrument

%ConsumptionVar: One of 'JW' (ND+Services) 'Hall' (only ND)

%RunRidge: true/false
%note: controls both zero-beta estimate and consumption
%note: if true, overrides GenSEs and won't generate standard errors

%GenSEs: true/false:
% note: controls both generation of SEs and production of S-test graph

%sigma_max: value for maximum sigma to test
%10 produces most graphs, 20 is sufficient to illustrate non-linear
%'rare disasters' effect

%SplitSample/SplitYear: fit regressions only on earlier data, plot full
%sample



opts = struct('Name','Main');
opts.AssetFile = '27_plus_Industry_Portfolios_Nominal.csv';
opts.VaryingBeta = false;
opts.LinearConsumption = false;
opts.NLConsFactor = true;
opts.SigmaInit = 5;
opts.Factors = {'Mkt','SMB','HML','RMW','CMA','term_spread','DEF'};
opts.Instruments = {'RF','UMP','EBP','TSP','CPI_rolling'};
opts.LaggedCons = false;
opts.ConsumptionVar = 'JW';
opts.RunRidge = false;
opts.GenSEs = true;
opts.sigma_max = 10;
opts.NoCOVID = false;
opts.har = 'COV';
opts.SplitSample = false;
optsDefault = opts;

%run tests expects variable named 'opts';
RunTest;


opts = optsDefault;
opts.Name = 'NoDrop20';
opts.AssetFile = 'NoDrop20/27_plus_Industry_Portfolios_Nominal.csv';
RunTest;

opts = optsDefault;
opts.Name = 'FF5Industry';
opts.AssetFile = 'FF5_plus_Industry_Portfolios_Nominal.csv';
RunTest;

opts = optsDefault;
opts.Name = 'Sigma1';
opts.SigmaInit = 1;
RunTest;

opts = optsDefault;
opts.Name = 'LinearCons';
opts.LinearConsumption = true;
opts.NLConsFactor = false;
RunTest;

opts = optsDefault;
opts.Name = 'NoCons';
opts.LinearConsumption = false;
opts.NLConsFactor = false;
RunTest;

opts = optsDefault;
opts.Name = 'MktOnly';
opts.Factors = {'Mkt'};
RunTest;

opts = optsDefault;
opts.Name = 'FF3Only';
opts.Factors = {'Mkt','SMB','HML'};
RunTest;

opts = optsDefault;
opts.Name = 'AltBAAS';
opts.Instruments = {'RF','UMP','TSP','CPI_rolling','BAAS'};
RunTest;

opts = optsDefault;
opts.Name = 'AltCAPE';
opts.Instruments = {'RF','UMP','EBP','CAPE','TSP','CPI_rolling'};
RunTest;

opts = optsDefault;
opts.Name = 'AltDP';
opts.Instruments = {'RF','UMP','EBP','DP_ratio','CAPE', 'TSP','CPI_rolling'};
RunTest;

opts = optsDefault;
opts.Name = 'LagCons';
opts.LaggedCons = true;
RunTest;

opts = optsDefault;
opts.Name = 'Shadow';
opts.Instruments = {'RF','UMP','EBP','TSP','CPI_rolling','shadow_spread'};
RunTest;

opts = optsDefault;
opts.Name = 'VaryingBetas';
opts.VaryingBeta = true;
RunTest;

opts = optsDefault;
opts.Name = 'NDOnly';
opts.ConsumptionVar = 'Hall';
RunTest;

opts = optsDefault;
opts.Name = 'Ridge';
opts.RunRidge = true;

RunTest;

opts = optsDefault;
opts.Name = 'NoCOVID';
opts.NoCOVID = true;
RunTest;

opts = optsDefault;
opts.Name = 'Pre2004';
opts.SplitSample = true;
opts.SplitYear = 2004;
opts.RunRidge = true;
RunTest;

opts = optsDefault;
opts.Name = 'Pre1989';
opts.SplitSample = true;
opts.SplitYear = 1989;
opts.RunRidge = true;
RunTest;



