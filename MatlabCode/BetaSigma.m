function [Beta, Sigma,alphas] = BetaSigma(R, Rm, Z, Rb, ConsG, inflation, iotaN, iotaM, Rf, gamma, sigma,NLConsFactor)
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
SigmaP = DirectKernel(eP');

%undo the pre-conditioning to get final cov matrix estimate
Sigma = SigmaFroot*SigmaP*SigmaFroot;

%code expects transposed shape
Beta = Beta';

