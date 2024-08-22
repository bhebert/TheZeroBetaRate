function [weight] = PortfolioWeight(Beta,Sigma,iotaN)

        %construct zero-beta portfolio
        %todo: figure out way to compute with backslash operator to speed up
        M = size(Beta,2);
        Amat = [iotaN, Beta];
        Sinv = inv(Sigma); %old: pinv
        e1 = [1; zeros(M,1)];
        mult = inv(Amat' * Sinv * Amat) * e1; %old: pinv
        weight = (Sinv* Amat * mult)';
  
end