function sigmahat = DirectKernel(X)
%DirectKernel - Computes the nonlinear shrinkage estimator of the covariance matrix.
%
% Syntax: sigmahat = DirectKernel(X)
%
% Input: X - N x P, N: number of observations, P: number of variables.
    

    % extract sample eigenvalues sorted in ascending order and eigenvectors
    [n,p] = size(X);
    sample = (X'*X)./n;
    [u, lambda] = eig(sample, 'vector');
    [lambda, isort] = sort(lambda);
    u = u(:, isort);

    % compute direct kernel estimator
    lambda = lambda(max(1 , p-n+1) : p);
    L = repmat(lambda,[1 min(p,n)]);
    h = n^(-0.35);
    ftilde=mean(sqrt(max(0,4*L'.^2*h^2-(L-L').^2))./(2*pi*L'.^2*h^2),2);
    Hftilde=mean((sign(L-L').*sqrt(max(0,(L-L').^2-4*L'.^2*h^2))-L+L')./(2*pi*L'.^2*h^2),2);

    if p <= n
        
        dtilde = lambda./((pi*(p/n)*lambda.*ftilde).^2 ...
            +(1-(p/n)-pi*(p/n)*lambda.*Hftilde).^2);

    else
        Hftilde0=(1-sqrt(1-4*h^2))/(2*pi*h^2)*mean(1./lambda);
        dtilde0=1/(pi*(p-n)/n*Hftilde0);
        dtilde1=lambda./(pi^2*lambda.^2.*(ftilde.^2+Hftilde.^2));
        dtilde=[dtilde0*ones(p-n,1);dtilde1];
    end

    dhat=pav(dtilde);
    sigmahat=u*(repmat(dhat,[1 p]).*u');
end