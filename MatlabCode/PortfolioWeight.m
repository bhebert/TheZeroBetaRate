function [weight] = PortfolioWeight(Beta,Sigma,iotaN)


    %construct zero-beta portfolio
    %todo: figure out way to compute with backslash operator to speed up
    M = size(Beta,1);
    Amat = [iotaN, Beta];
    Sinv = pinv(Sigma);
    e1 = [1; zeros(M,1)];
    mult = pinv(Amat' * Sinv * Amat) * e1;
    weight = (Sinv* Amat * mult)';


end