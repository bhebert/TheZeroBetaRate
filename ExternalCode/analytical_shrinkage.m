function [sigmatilde,dtilde]=analytical_shrinkage(X,k)

% this code is published in the supplemental materials of
% Ledoit, Olivier, and Michael Wolf. "Analytical nonlinear shrinkage of large-dimensional covariance matrices." 
% Ann. Statist. 48(5): 3043-3065 (October 2020). DOI: 10.1214/19-AOS1921

%the code was downloaded from Michael Wolf's website,
%https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html
%accessed Mar 8, 2023

% X is the raw data matrix of size n x p:
% - the rows correspond to observations
% - the columns correspond to variables

% If the second (optional) parameter k is absent, not-a-number, or empty, 
% the algorithm demeans the data by default, and then adjusts 
% the effective sample size accordingly by subtracting one.

% If the user inputs k = 0, then no demeaning takes place and
% the effective sample size remains n.

% If the user inputs k >= 1, then it signifies that the data X 
% has already been demeaned or otherwise pre-processed; for example,
% the data might constitute OLS residuals based on a linear regression
% model with k regressors. No further demeaning takes place then,
% but the effective sample size is adjusted accordingly by subtracting k.

[n,p]=size(X); % important: sample size n must be >= 12
if (nargin<2)||isnan(k)||isempty(k) % default setting
   X=X-repmat(mean(X),[n 1]); % demean the raw data matrix
   k=1; % subtract one degree of freedom 
end
n=n-k; % effective sample size
% extract sample eigenvalues sorted in ascending order and eigenvectors
sample=(X'*X)./n;
[u,lambda]=eig(sample,'vector');
[lambda,isort]=sort(lambda);
u=u(:,isort);
% compute analytical nonlinear shrinkage kernel formula
lambda=lambda(max(1,p-n+1):p);
L=repmat(lambda,[1 min(p,n)]);
h=n^(-1/3);
H=h*L';
x=(L-L')./H;
ftilde=(3/4/sqrt(5))*mean(max(1-x.^2./5,0)./H,2);
Hftemp=(-3/10/pi)*x+(3/4/sqrt(5)/pi)*(1-x.^2./5) ...
   .*log(abs((sqrt(5)-x)./(sqrt(5)+x)));
Hftemp(abs(x)==sqrt(5))=(-3/10/pi)*x(abs(x)==sqrt(5));
Hftilde=mean(Hftemp./H,2);
if p<=n
   dtilde=lambda./((pi*(p/n)*lambda.*ftilde).^2 ...
      +(1-(p/n)-pi*(p/n)*lambda.*Hftilde).^2);
else
   Hftilde0=(1/pi)*(3/10/h^2+3/4/sqrt(5)/h*(1-1/5/h^2) ...
      *log((1+sqrt(5)*h)/(1-sqrt(5)*h)))*mean(1./lambda); dtilde0=1/(pi*(p-n)/n*Hftilde0);
   dtilde1=lambda./(pi^2*lambda.^2.*(ftilde.^2+Hftilde.^2)); 
   dtilde=[dtilde0*ones(p-n,1);dtilde1];
end
% compute analytical nonlinear shrinkage estimator
sigmatilde=u*diag(dtilde)*u';