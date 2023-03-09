function momjacbobian = InstMomentsFullJacobian(Theta,Beta,alphas, R, Rm, Z, Rb, iotaN, iotaM, ConsG,inflation, weight, sig,NLConsFactor,varargin)

%  Computes the full jacobian of all moments
%  Takes a mixed approach.
%  For projections moments, the jacobian is computed analytically.
%  For asset pricing and euler moments, the jacobian is computed numerically.

    K = size(Z, 1);
    

    jacProj = InstMomentsProjectionJacobian(Theta,Beta,alphas, R, Rm, Z, Rb, iotaN, iotaM, ConsG,inflation,NLConsFactor);
    
    xstart = [Theta(1:end-1);reshape(Beta,size(Beta,1)*size(Beta,2),1);alphas'];
    

    %reference for what function call should look like
    %use 'ZB' and not 'ZBFull' because we don't want proj moments
    %okay if Rfex = [] because it won't be used
    %InstMomentsConc([Theta1;sig],alphas,Beta,weight, Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, ConsG,inflation,Rfex,'ZB',NLConsFactor);


    func = @(x) InstMomentsConc([x(1:2+K);sig], ...
            reshape(x(3+K+size(Beta,1)*size(Beta,2):end),size(Beta,1),1), ...
            reshape(x(3+K:2+K+size(Beta,1)*size(Beta,2)), ...
            size(Beta,1),size(Beta,2)), weight, ...
            R,Rm,Z,Rb,iotaN,iotaM,ConsG,inflation,[],'ZB',NLConsFactor);

    
    jacAP = jacobianest(func, xstart);

    momjacbobian = [jacProj; jacAP];

end

