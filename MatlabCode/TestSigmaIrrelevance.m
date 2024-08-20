%test script to verify irrelevance of sigma value
%this is essentially a test of the numerical method

clear;
load("../Output/MainRun_Main.mat");

[~, Theta2] = InstrumentGMMWrapperConc(Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, cons_gr_ann/12, inflation, Rfex, psi,2,'ZB',opts.har,opts.NLConsFactor,opts.SigmaType);
[~, Theta10] = InstrumentGMMWrapperConc(Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, cons_gr_ann/12, inflation, Rfex, psi,10,'ZB',opts.har,opts.NLConsFactor,opts.SigmaType);

disp("the two vectors should be equal except for the last two params")
disp("which are the discount rate and the sigma value")
Theta2 - Theta10