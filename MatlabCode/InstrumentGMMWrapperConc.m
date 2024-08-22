function [test, Theta1_final, Sv, portRet, weight, Sigma, Beta, alphas, mmoments,cmoments,amoments,thresh,mvar,pvar] = InstrumentGMMWrapperConc(Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, ConsG,inflation,Rfex,psi,sig,asset,har,NLConsFactor,SigmaType)
% This function is a wrapper for the GMM optimization process.



%-----------------------------------------------------------------
%                           Estimation
%-----------------------------------------------------------------

N = size(Rinput, 1); T = size(Rinput, 2); K = size(Zinput, 1);


%a better guess for rho_init
rho_init = (0.08 - sig*0.015)/12;

if strcmp(asset,'RF') || strcmp(asset,'Mkt')
    Theta0_init = [rho_init];
    ub = 0.02;
    lb = -ub;
    TypX = 0.01;
    %don't do global search if failing for some reason
    cnt = 7;
else
   
    % run once to build initial guess with zbrate = safe rate
    [Beta, Sigma] = BetaSigma(Rinput, Rminput, Zinput, Rbinput, ConsG, inflation,iotaN, iotaM, 0, zeros(K,1), sig,NLConsFactor,SigmaType);
    weight =  PortfolioWeight(Beta,Sigma,iotaN);
    portRet_init = weight * Rinput; 

    %use infeasible OLS as starting point
    Theta0 = [ones(1,size(Zinput,2)); Zinput]' \ (portRet_init-Rbinput)';

    %repeat to iterate on a better guess
    [Beta, Sigma] = BetaSigma(Rinput, Rminput, Zinput, Rbinput, ConsG, inflation,iotaN, iotaM, Theta0(1),reshape(Theta0(2:1+K), K, 1), sig,NLConsFactor,SigmaType);
    weight =  PortfolioWeight(Beta,Sigma,iotaN);
    portRet_init = weight * Rinput; 

    Theta0 = [ones(1,size(Zinput,2)); Zinput]' \ (portRet_init-Rbinput)';
    Theta0_init = [Theta0; rho_init];

    %bounds at large values to avoid singular matrices
    ub = [ones(K+1,1)*3*max(abs(Theta0_init));.02];
    lb = -ub;
    TypX = [ones(K+1,1);0.01];
    cnt=1;
end

%options_min = optimoptions('fminunc', 'Display','iter', 'MaxIterations', 1000,  'MaxFunctionEvaluations', 1e6,'OptimalityTolerance', 10^-8,'StepTolerance',10^-14,'FunctionTolerance',10^-10);
options_min = optimoptions('fmincon', 'Display','off', 'MaxIterations', 100000,  'MaxFunctionEvaluations', 1e10,'OptimalityTolerance', 10^-10,'StepTolerance',10^-14,'FunctionTolerance',10^-10,'TypicalX',TypX);

% GMM (psi>0: with ridge penalty)
f1 = @(Theta) InstrumentGMMConcOpt([Theta;sig], Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, ConsG,inflation,Rfex,asset,NLConsFactor,SigmaType) + psi * sum(Theta(2:end-2).^2); 

%[Theta1, ~]  = fminunc(f1, Theta0_init, options_min);
[Theta1, obj, flag]  = fmincon(f1, Theta0_init, [],[],[],[], lb,ub,[], options_min);

%objective will not be zero if running ridge, so can't use objective to
%detect error
while (psi == 0 && obj > 10^-6 && cnt < 6)
    %try again
    %use fixed seed
    rng(12302018,'twister');
    gs = GlobalSearch('FunctionTolerance',10^-10); %MultiStart('FunctionTolerance',10^-8,'UseParallel',true);

    %iterate 10 times to move away from the stuck point
    for i=1:10
        [Beta, Sigma] = BetaSigma(Rinput, Rminput, Zinput, Rbinput, ConsG, inflation,iotaN, iotaM, Theta1(1), reshape(Theta1(2:1+K), K, 1), sig,NLConsFactor,SigmaType);
        weight =  PortfolioWeight(Beta,Sigma,iotaN);
        portRet_init = weight * Rinput; 

        %use infeasible OLS as starting point
        Theta0_1 = [ones(1,size(Zinput,2)); Zinput]' \ (portRet_init-Rbinput)';
        Theta1 = [Theta0_1;Theta1(end)];
    end

    %try a relatively narrow region near this guess
    ub = Theta1 + [0.02*ones(K+1,1);.001];
    lb = Theta1 - [0.02*ones(K+1,1);.001];

    problem = createOptimProblem('fmincon','x0',Theta1,'objective',f1,'lb',lb,'ub',ub); %'options',options_min
    %'lb',lb,'ub',ub
    
    [Theta1, obj]  = run(gs,problem); % fmincon(f1, Theta0_init, [],[],[],[], -inf*ones(size(Theta0_init)),inf*ones(size(Theta0_init)),[], options_min2);
    cnt = cnt + 1;
    
end

if (psi==0 && obj>10^-6) || flag < 1

    obj
    Theta1
    error('Numerical algorithm did not converge');
end



Theta1_final = [Theta1;sig];

if strcmp(asset,'RF') || strcmp(asset,'Mkt')
    Beta = [];
    Sigma = [];
    alphas = [];
    weight = [];
else
    %the first parameter is the constant
    Rfcons = Theta1(1);
    gamma = reshape(Theta1(2:1+K), K, 1);

    %run to get some intermediate results
    [Beta, Sigma, alphas] = BetaSigma(Rinput, Rminput, Zinput, Rbinput, ConsG, inflation,iotaN, iotaM, Rfcons, gamma, sig,NLConsFactor,SigmaType);
    weight =  PortfolioWeight(Beta,Sigma,iotaN);
end

[mmoments,amoments, cmoments, ~,portRet] = InstMomentsConc([Theta1;sig],alphas,Beta,weight, Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, ConsG,inflation,Rfex,asset,NLConsFactor);



Nm = size(amoments,1)+1;

%don't make standard errors if running ridge
if psi == 0

    % the jacobian with respect to the discount rate can be explicitly 
    % computed for these cases
    if strcmp(asset,'RF') || strcmp(asset,'Mkt')
        
        Om = Wopt(cmoments,har)/T;
        e1 = [1;zeros(K,1)];

        % can analytically compute jacobian in this case
        % the moment errors are scaled by 100 for numerical reasons
        jacRF = -100*(e1+mmoments(Nm:end)/100);

        %Wmat in this case is e1 * e1'
        % general formula: GWGiGW = (jac'*Wmat*jac) \ (jac' * Wmat);   
        % in this case, 
        GWGiGW = e1' / jacRF(1);
        
        %the usual variance formula for parameters, not actually used in
        %anything
        pvar = GWGiGW * Om * (GWGiGW)';

        %the usual GMM moment variance formula
        %the first moment has no variance due to exact identification, omit
        mvar = (eye(size(Om)) - jacRF*GWGiGW) * Om * (eye(size(Om)) - jacRF*GWGiGW)'; 
        mvar = mvar(2:end,2:end);
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




