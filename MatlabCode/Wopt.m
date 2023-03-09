function [omega_hat_mat] = Wopt(gs,type)
% variance matrix estimator
% in main specification, use no lags
% per Hansen Singleton 1983, Cochrane 2005 chapter 11.7

% As an alternative, use Newey-West 1987
% NW code from 
% Lazarus, E., Lewis, D. J., Stock, J. H., & Watson, M. W. (2018). HAR inference: Recommendations for practice. Journal of Business & Economic Statistics, 36(4), 541-559.
% the NW currently uses the short (textbook) NW tuning
%   code adapted from their file "har_ols.m" and "hac_ols.m"

T = size(gs,2);

%center for varcov esimtation ala Hall 2000
gsm = gs - mean(gs,2,'omitnan');


if contains(type,'COV')
    omega_hat_mat = cov(gs','partialrows');

else

    nma = ceil(0.75*T^(1/3));


    omega_hat_mat=zeros(size(gs,1),size(gs,1));
    z = gsm';
   

    % Form Kernel 

    kern = 1 - (0:1:nma)/(nma+1);


    % Form Hetero-Serial Correlation Robust Covariance Matrix 
    for ii = -nma:nma
        if ii <= 0 
            r1=1; 
            r2=size(z,1)+ii; 
        else 
            r1=1+ii; 
            r2=size(z,1); 
        end

        omega_hat_mat=omega_hat_mat + kern(abs(ii)+1)*(z(r1:r2,:)'*z(r1-ii:r2-ii,:))/T;
    end

end


end