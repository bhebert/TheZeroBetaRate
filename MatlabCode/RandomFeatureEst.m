function [Theta,portRet3,p_cons,norm,weight,cbetas] = RandomFeatureEst(Rinput, Rminput, Zinput, Rbinput,cons_gr,inflation,iotaN, iotaM, sigma,NLConsFactor,psi,psic)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    thresh = 10^-8;
    Rf = 0;
    gamma =  zeros(size(Zinput,1),1);
    i=0;
    norm = thresh*10;

    Zmat = Zinput * Zinput';
    Imat = inv(Zmat+psi*eye(size(Zmat)));
    Imatc = inv(Zmat+psic*eye(size(Zmat)));

    while i < 100 && norm > thresh
    
          % run once to build initial guess with zbrate = safe rate
        [Beta, Sigma] = BetaSigma(Rinput, Rminput, Zinput, Rbinput, cons_gr, inflation,iotaN, iotaM, Rf, gamma, sigma,NLConsFactor);
        weight =  PortfolioWeight(Beta,Sigma,iotaN);
        portRet3 = weight * Rinput; 

        %use infeasible OLS as starting point
        %Theta0 = ridge(portRet3'-Rbinput',Zinput',psi,0);
        gamma = Imat * Zinput * (portRet3'-Rbinput' - mean(portRet3'-Rbinput'));
        Rf = mean(portRet3'-Rbinput') - mean(gamma'*Zinput);
        zbrate = gamma'*Zinput + Rf + Rbinput; % zero-beta rate
        beta = (Rminput-iotaM*zbrate)' \ (portRet3-zbrate)';
        
        norm = sqrt(beta' * beta);
        i = i+1;
    end
  
    
    cbetas = Imatc * Zinput * (cons_gr-mean(cons_gr));
    cmean = mean(cons_gr) - mean(cbetas'*Zinput);
    p_cons = 12*(cmean + cbetas'*Zinput);
    cbetas = [cmean; cbetas];
    Theta = [Rf;gamma; -.02/12; sigma];

end