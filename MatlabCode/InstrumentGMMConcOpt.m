function [obj] = InstrumentGMMConcOpt(Theta, R, Rm, Z, Rb, iotaN, iotaM, ConsG,inflation,Rfex,asset,NLConsFactor)
%wrapper function to compute GMM objective
T = size(R, 2);

K = size(Z, 1);
Rfcons = Theta(1);
gamma = reshape(Theta(2:1+K), K, 1);
sigma = Theta(end);

if strcmp(asset,'ZB')
    [Beta, Sigma] = BetaSigma(R, Rm, Z, Rb, ConsG, inflation,iotaN, iotaM, Rfcons, gamma, sigma,NLConsFactor);
    weight =  PortfolioWeight(Beta,Sigma,iotaN);
    Nm = size(Z,1) + 1;
else
    Beta = [];
    weight = [];
    Nm = 1;
end

mmoments = InstMomentsConc(Theta, Beta, weight, R, Rm, Z, Rb, iotaN, iotaM, ConsG,inflation,Rfex, asset,NLConsFactor);

obj = T*mmoments(1:Nm)'* mmoments(1:Nm);


