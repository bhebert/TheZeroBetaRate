function [Beta, Sigma,alphas] = BetaSigma(R, Rm, Z, Rb, ConsG, inflation, iotaN, iotaM, Rf, gamma, sigma,NLConsFactor,SigmaType)
% This function recover beta via OLS and recover Sigma via Direct Kernel
% Method.


T = size(R, 2);

% ---- OLS to Recover Beta ---- %
zbrate = gamma'*Z + Rf + Rb; % zero-beta rate

X = R - iotaN*zbrate; % excess returns of the portfolios

Xm = Rm - iotaM*zbrate; % excess returns of the factors

if NLConsFactor
    % consumption factor to orthogonalize against
    ConsF = 100*(exp(-sigma*ConsG/100 - inflation')-1);

    Xv = [Xm;ConsF';ones(1,T)]';
    %covariance matrix of factors
    SigmaM = cov([Xm;ConsF']');
else
    Xv = [Xm;ones(1,T)]';
    
    SigmaM = cov(Xm');
end

BetaFull = Xv \ X';

Beta = BetaFull(1:end-1,1:end);
alphas = BetaFull(end,1:end);

eF = X' - Xv * BetaFull;


%covariance matrix of residuals
SigmaE = cov(eF);

if strcmp(SigmaType,'Sample')
    %test code to use standard covariance estimate in place of LW
    Sigma = Beta'*SigmaM*Beta + SigmaE;
elseif strcmp(SigmaType,'I')
    Sigma = eye(size(SigmaE));
elseif strcmp(SigmaType,'Diag')
    Sigma = Beta'*SigmaM*Beta + diag(diag(SigmaE));
elseif strcmp(SigmaType,'PCA')
    [coeff,latent] = pcacov(SigmaE);
    SigmaE2 = coeff(:,1:3)*diag(latent(1:3))*coeff(:,1:3)';
    Sigma = Beta'*SigmaM*Beta + SigmaE2 + diag(diag(SigmaE-SigmaE2));
else

    %covariance matrix implied by factor model
    %use this to precondition before LW2017 estimator
    %see LW2017 section 4.2
    SigmaF = Beta'*SigmaM*Beta + diag(diag(SigmaE));
    
    %first take symmteric square root
    SigmaFroot = sqrtm(SigmaF);
    
    %then de-mean excess returns and multiply by
    %inverse of square root
    mu = mean(X,2);
    eP = SigmaFroot\(X-mu);
    
    %use LW2017 estimator with direct kernel method to get
    %regularized covariance matrix
    %SigmaP = DirectKernel(eP');
    
    %use LW2020 estimator with analytical method to get
    %regularized covariance matrix
    SigmaP = analytical_shrinkage(eP');
    
    %undo the pre-conditioning to get final cov matrix estimate
    Sigma = SigmaFroot*SigmaP*SigmaFroot;
end


%code expects transposed shape
Beta = Beta';

