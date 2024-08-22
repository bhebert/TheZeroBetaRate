%test script to verify irrelevance of sigma value
%this is essentially a test of the numerical method

clear;
load("../Output/MainRun_MktOnly.mat");
opts.SigmaType="LW";
addpath("../ExternalCode");

warning('error','MATLAB:nearlySingularMatrix');

[~, Theta2,~, portRet2,] = InstrumentGMMWrapperConc(Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, cons_gr_ann/12, inflation, Rfex, psi,2,'ZB',opts.har,opts.NLConsFactor,opts.SigmaType);
[~, Theta10,~,portRet10] = InstrumentGMMWrapperConc(Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, cons_gr_ann/12, inflation, Rfex, psi,10,'ZB',opts.har,opts.NLConsFactor,opts.SigmaType);

disp("the two vectors should be equal except for the last two params")
disp("which are the discount rate and the sigma value")
Theta2 - Theta10

Rf2 = Theta2(1);
gamma2 = reshape(Theta2(2:1+K), K, 1);

Rf10 = Theta10(1);
gamma10 = reshape(Theta10(2:1+K), K, 1);

zbrate2 = gamma2'*Zinput + Rf2 + Rbinput; % zero-beta rate
zbrate10 = gamma10'*Zinput + Rf10 + Rbinput; % zero-beta rate

mean(portRet2-zbrate2)

(portRet2-zbrate2)*(Rminput-iotaM*zbrate2)'/T

([Zinput;ones(size(portRet3))]' \ (portRet2-Rbinput)')' - [gamma2;Rf2]'



mean(portRet10-zbrate10)

(portRet10-zbrate10)*(Rminput-iotaM*zbrate10)'/T

([Zinput;ones(size(portRet3))]' \ (portRet10-Rbinput)')' - [gamma10;Rf10]'
