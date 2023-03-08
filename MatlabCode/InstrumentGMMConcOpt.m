function [obj] = InstrumentGMMConcOpt(Theta, R, Rm, Z, Rb, iotaN, iotaM, ConsG,inflation,Rfex,asset,NLConsFactor)
%wrapper function to compute GMM objective
T = size(R, 2);

if strcmp(asset,'ZB')
    [Beta, Sigma] = BetaSigma(R, Rm, Z, Rb, ConsG, inflation,iotaN, iotaM, 0, zeros(K,1), sigma,NLConsFactor);
    weight =  PortfolioWeight(Beta,Sigma,iotaN);
    Nm = size(Z,1) + 1;
else
    Beta = [];
    weight = [];
    Nm = 1;
end

mmoments = InstMomentsConc(Theta, Beta, weight, R, Rm, Z, Rb, iotaN, iotaM, ConsG,inflation,Rfex, asset,NLConsFactor);

obj = T*mmoments(1:Nm)'* mmoments(1:Nm);


