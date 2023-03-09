function [mmoments, amoments, consmoments, projmoments, portRet] = InstMomentsConc(Theta, alphas, Beta, weight, R, Rm, Z, Rb, iotaN, iotaM, ConsG,inflation,Rfex, asset,NLConsFactor)

N = size(R, 1); T = size(R, 2);  K = size(Z, 1);

if length(Theta) > 2
    Rf = Theta(1);
    gamma = reshape(Theta(2:1+K), K, 1);
end

sigma = Theta(end);

rho = Theta(end-1);


if NLConsFactor
    M = size(Rm, 1)+1;
else
    M = size(Rm,1);
end

consSDF = exp(-(rho + sigma*ConsG'/100 +inflation));
zmat = [ones(1,T);Z];

if strcmp(asset,'RF')

    portRet = Rfex;
    amoments = [];
    projmoments = [];

    %Euler equation error
    conserr = 100*(consSDF.*(1+portRet/100)-1);


elseif strcmp(asset,'Mkt')
    portRet = Rm(1,:);
    amoments = [];
    projmoments = [];

    %Euler equation error
    conserr = 100*(consSDF.*(1+portRet/100)-1);
    

else
    zbrate = gamma'*Z + Rf + Rb; % zero-beta rate


    %orthogonal projection matrix
    Hmat = eye(N) - Beta*pinv(Beta);

    %note that weight*Hmat = weight

    %compute return of portfolio
    portRet = weight * Hmat * R; 


    conserr = 100*(consSDF.*(1+zbrate/100)-1);


    

    %moments are product of surprise and instruments
    amoments = (weight * Hmat * (R - zbrate)).*zmat;

    %don't compute this unless computing standard errors
    if strcmp(asset,'ZBFull')
        if NLConsFactor 
            ConsF = 100*(exp(-sigma*ConsG/100 - inflation')-1);
            e = R - alphas' - (iotaN - Beta*[iotaM;0])*zbrate - Beta*[Rm;ConsF'];
            Fmat = [ones(1,T);Rm-iotaM*zbrate;ConsF'];
        else
            e = R - alphas' - (iotaN - Beta*iotaM)*zbrate - Beta*Rm;
            Fmat = [ones(1,T);Rm-iotaM*zbrate];
        end

        e2 = permute(repmat(e,[1,1,M+1]),[2,1,3]);
        Fmat2 = permute(repmat(Fmat,1,1,N),[2,3,1]);

        projmoments = reshape(e2 .* Fmat2,T,[]);
    else
        projmoments = [];
    end

end



%moments are product of Euler error and instruments
consmoments = repmat(conserr,K+1,1) .* zmat;
projmoments = projmoments';
momfull=[projmoments; amoments; consmoments];

mmoments = mean(momfull,2);

end