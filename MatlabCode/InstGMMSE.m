function [pvar,mvar] = InstGMMSE(Theta,Beta,alphas, R, Rm, Z, Rb, iotaN, iotaM, ConsG,inflation,weight,har,sig,NLConsFactor)

    K = size(Z, 1);
    T= size(Z,2);

 
    %don't need Rfex in this call
    [~,amoments, cmoments, projmoments] = InstMomentsConc(Theta,alphas,Beta,weight, R, Rm, Z, Rb, iotaN, iotaM, ConsG,inflation,[],'ZBFull',NLConsFactor);
    momfull=[projmoments; amoments; cmoments];

    %[~, momfull] = InstMomentsFull(Theta,Beta,alphas', R, Rm, Z, Rb, iotaN, iotaM, ConsG,inflation,weight,NLConsFactor);
    

   

    Nm = size(momfull,1);

    %this is the exact identifiction weight matrix
    Wmat = zeros(Nm);
    Wmat(1:(Nm-K),1:(Nm-K)) = eye(Nm-K);


    Omega = Wopt(momfull,har);


    %the below commented code would compute the jacobian numerically
    %we have instead written a function to compute the projection moment
    %part of the jacobian analytically, which speeds things up
    %we will compute the asset pricing/euler part numerically
    %xinit = [Theta(1:end-1);reshape(Beta,size(Beta,1)*size(Beta,2),1);alphas'];

%     func = @(x) InstMomentsFull([x(1:2+K);sig],reshape(x(3+K:2+K+size(Beta,1)*size(Beta,2)), ...
%         size(Beta,1),size(Beta,2)),reshape(x(3+K+size(Beta,1)*size(Beta,2):end),size(Beta,1),1), ...
%         R,Rm,Z,Rb,iotaN,iotaM,ConsG,inflation,weight);
% 
%     jac = jacobianest(func,xinit);
    

    jac = InstMomentsFullJacobian(Theta,Beta,alphas, R, Rm, Z, Rb, iotaN, iotaM, ConsG,inflation, weight, sig,NLConsFactor);



    %the usual GMM parameter variance formula
    GWGiGW = (jac'*Wmat*jac) \ (jac' * Wmat);   
    pvar = 1/T* GWGiGW * Omega * (GWGiGW)';

    %the usual GMM moment variance formula
    mvar = 1/T * (eye(Nm) - jac*GWGiGW) * Omega * (eye(Nm) - jac*GWGiGW)'; 

    %these are the only non-zero parts
    mvar = mvar(Nm-K+1:end,Nm-K+1:end);

    %only care about uncertainty with respect to non-alpha/beta parameters
    pvar = pvar(1:K+2,1:K+2);

end