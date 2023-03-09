function [test, Theta1_final, Sv, portRet, weight, Sigma, Beta, alphas, mmoments,cmoments,amoments,thresh,mvar,pvar] = InstrumentGMMWrapperConc(Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, ConsG,inflation,Rfex,psi,sig,asset,har,NLConsFactor)
% This function is a wrapper for the GMM optimization process.



%-----------------------------------------------------------------
%                           Estimation
%-----------------------------------------------------------------

N = size(Rinput, 1); T = size(Rinput, 2); K = size(Zinput, 1);

%an initial discount rate of about 2%/yr
rho_init = -.02/12;

if strcmp(asset,'RF') || strcmp(asset,'Mkt')
    Theta0_init = [rho_init];
else
   
    % run once to build initial guess with zbrate = safe rate
    [Beta, Sigma] = BetaSigma(Rinput, Rminput, Zinput, Rbinput, ConsG, inflation,iotaN, iotaM, 0, zeros(K,1), sig,NLConsFactor);
    weight =  PortfolioWeight(Beta,Sigma,iotaN);
    portRet_init = weight * Rinput; 

    %use infeasible OLS as starting point
    Theta0 = [ones(1,size(Zinput,2)); Zinput]' \ (portRet_init-Rbinput)';
    Theta0_init = [Theta0; rho_init];
end


options_min = optimoptions('fminunc', 'Display','off', 'MaxIterations', 1000,  'MaxFunctionEvaluations', 1e6,'OptimalityTolerance', 10^-8);




% GMM (psi>0: with ridge penalty)
f1 = @(Theta) InstrumentGMMConcOpt([Theta;sig], Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, ConsG,inflation,Rfex,asset,NLConsFactor) + psi * sum(Theta(2:end-2).^2); 

[Theta1, ~]  = fminunc(f1, Theta0_init, options_min);

Theta1_final = [Theta1;sig];

%the first parameter is the constant
Rfcons = Theta(1);
gamma = reshape(Theta(2:1+K), K, 1);

%run to get some intermediate results
[Beta, Sigma, alphas] = BetaSigma(Rinput, Rminput, Zinput, Rbinput, ConsG, inflation,iotaN, iotaM, Rfcons, gamma, sig,NLConsFactor);
weight =  PortfolioWeight(Beta,Sigma,iotaN);
[mmoments,amoments, cmoments, ~,portRet] = InstMomentsConc([Theta1;sig],Beta,weight, Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, ConsG,inflation,Rfex,asset,NLConsFactor);



Nm = size(amoments,1)+1;

%don't make standard errors if running ridge
if psi == 0

    % the jacobian with respect to the discount rate can be explicitly 
    % computed for these cases
    if strcmp(asset,'RF') || strcmp(asset,'Mkt')
        Om = Wopt(cmoments,har)/T;
        e1 = [1;zeros(K,1)];
        mat = eye(size(Om)) - e1 * (e1+mmoments(Nm:end))';
        mvar = mat' * Om * mat;
        mvar = mvar(2:end,2:end);
        pvar = Om(1,1);
    else

        [pvar,mvar] = InstGMMSE([Theta1;sig],Beta,alphas, Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, ConsG,inflation,weight,har,sig,NLConsFactor);

    end
end

if psi == 0
    %compute S-statistic
    Sv = mmoments(Nm+1:end)' * (mvar \ mmoments(Nm+1:end));

    %95% threshold
    thresh = chi2inv(0.95,K);
    test = Sv < thresh;
else
    pvar = [];
    mvar = [];

    Sv = 0;
    thresh = 0;
    test = 0;

end



end




