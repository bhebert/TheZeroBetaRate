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

%NoCOVID: true/false, end data in Dec 2019

%SplitSample/SplitYear: fit regressions only on earlier data, plot full
%sample



opts = struct('Name','Main');

%June 2023: switched to also include investment/profitability 3x3x3 sorts
%these sorts now make sense given fix to portfolio construction
%which eliminated problem of empty or almost-empty buckets
opts.AssetFile = 'FF5_plus_Industry_Portfolios_Nominal.csv';
opts.FactorFile = 'Factors_Nominal.csv';

opts.VaryingBeta = false;
opts.LinearConsumption = false;

%June 2023: switched main specification to not include cons factor
%makes exposition clearer
opts.NLConsFactor = false;

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
opts.RandomFeatures = false;
optsDefault = opts;

%run tests expects variable named 'opts';
RunTest;


opts = optsDefault;
opts.Name = 'NoDrop20';
opts.AssetFile = 'NoDrop20/FF5_plus_Industry_Portfolios_Nominal.csv';
RunTest;

%robustness test of portfolios constructed from BLLP data
%not included in main paper
%github code will not generate this file, so commented by default
%opts = optsDefault;
%opts.Name = 'AltDataBLLP';
%opts.AssetFile = 'AltData.csv';
%RunTest;


%June 2023: switched main spec to FF5, so robustness is not FF3
opts = optsDefault;
%opts.Name = 'FF5Industry';
%opts.AssetFile = 'FF5_plus_Industry_Portfolios_Nominal.csv';
opts.Name = 'FF3Industry';
opts.AssetFile = '27_plus_Industry_Portfolios_Nominal.csv';
RunTest;

%alternative portfolios by RA attempting to replicate 
%not included in paper, so commented by default
%github code will not generate the required input file
%opts = optsDefault;
%opts.Name = 'FF3IndustryPaul';
%opts.AssetFile = '27_plus_Industry_Portfolios_Nominal_paul.csv';
%RunTest;

%June 2023: Renamed this from Sigma1 to WithConsSigma1 due to switch in main spec 
opts = optsDefault;
opts.Name = 'WithConsSigma1';
opts.SigmaInit = 1;
opts.LinearConsumption = false;
opts.NLConsFactor = true;
RunTest;

opts = optsDefault;
opts.Name = 'LinearCons';
opts.LinearConsumption = true;
opts.NLConsFactor = false;
RunTest;

%June 2023: Renamed this from NoCons to WithCons due to switch in main spec 
opts = optsDefault;
opts.Name = 'WithConsSigma5';
opts.SigmaInit = 5;
opts.LinearConsumption = false;
opts.NLConsFactor = true;
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

%only run this if you have the informative factors data from
% Amman, Hemauer, and Straumann, "Informative Value, Profitability, and
% Investment Factors"
% Note that the file from those authors neededs to be processed before use.
% 1. Cut pre-Jan1973 data
% 2. Convert date format to same as our factors file
% 3. Add in term_spread and DEF from our file
% 4. Rename MP (market) to Mkt
% 5. Multiply the Mkt by 100
% 6. Add back in the t-bill rate (RF from the instrument file, one lag)
%    to the market-- the Mkt should be raw, not excess, return
% 6a. alternatively, use our Mkt variable-- after 5+6 the two are very
% close
% 7. Save as .csv

% opts = optsDefault;
% opts.Name = 'InformativeFactors';
% opts.Factors = {'Mkt','SMB','HML','RMW','CMA','term_spread','DEF'};
% opts.NoCOVID = true;
% opts.FactorFile = 'InformativeFactors.csv';
% RunTest;

%new random features procedure
opts = optsDefault;
opts.Name = 'RidgeRandom';
opts.RunRidge = true;
opts.RandomFeatures=true;
opts.RandomC = 10;
opts.RandomPsis = 100*[0.25 .5 1 1.5 2 2.5 3];
opts.RandomPsics = 10*[5 10 15 18 20 25 30];
opts.Randomgms = [0.25 0.5 1 2 4];
%opts.Randomgms = [0.5 0.6 0.7 0.8 0.9 1];
%opts.Randomgms = [1 1];
RunTest;


% Some split sample options. Did not make the final cut for the paper 
% but should still run if uncommented.
% opts = optsDefault;
% opts.Name = 'Pre2004';
% opts.SplitSample = true;
% opts.SplitYear = 2004;
% opts.RunRidge = true;
% RunTest;
% 
% opts = optsDefault;
% opts.Name = 'Pre1989';
% opts.SplitSample = true;
% opts.SplitYear = 1989;
% opts.RunRidge = true;
% RunTest;



