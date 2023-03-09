function momjacobian = InstMomentsProjectionJacobian(Theta,Beta,alphas, R, Rm, Z, Rb, iotaN, iotaM, ConsG,inflation,NLConsFactor)
%%
% 
%  Compute the Jacobian Matrix with respect to the projection moments.
%  
% 


    N = size(R, 1); T = size(R, 2);  K = size(Z, 1);

    if length(Theta) > 2
        Rf = Theta(1);
        gamma = reshape(Theta(2:end-2), K, 1);
    end

    sigma = Theta(end);
    rho = Theta(end-1);


    ConsF = 100*(exp(-sigma*ConsG/100 - inflation')-1);
    zbrate = gamma'*Z + Rf + Rb; % zero-beta rate
    
    % generate new factors and associated iotaM variables
    % just to make the expression for concise.
    if NLConsFactor
        iotaMp = [iotaM;0];
        Rmp = [Rm;ConsF'];
        M = size(Rm, 1)+1;
    else
        iotaMp = iotaM;
        Rmp = Rm;
        M = size(Rm, 1);
    end

    e = R - alphas' - (iotaN - Beta*iotaMp)*zbrate - Beta*Rmp;

    % dg(e)/dRf constant moment
    dge_dRf = zeros(N, 1);
    
    for i = 1:N
        dge_dRf(i) = - iotaN(i) + Beta(i, :)*iotaMp;
    end

    % dg(e)/dgamma constant moment
    dge_dgamma = zeros(N, K);

    for i = 1:N
        for j = 1:K
            dge_dgamma(i,j) = sum((- iotaN(i) + Beta(i, :)*iotaMp)*Z(j, :));
        end
    end
    dge_dgamma = dge_dgamma/T;

    % dg(e)/drho constant moment
    dge_drho = zeros(N, 1);


    % dg(e)/dB Constant moment
    dge_dB = zeros(N, N*M);
    
    for i = 1:N
        for p = 1:N
            for q = 1:M
                if i == p
                    dge_dB(i, N*(q-1)+p) = sum(iotaMp(q)*zbrate - Rmp(q,:));
                end
            end
        end
    end

    dge_dB = dge_dB/T;


    % dg(e)/dalpha Constant moment
    dge_dalpha = zeros(N, N);

    for i = 1:N
        for j = 1:N
            if i==j
                dge_dalpha(i,j) = -1;
            end
        end
    end

    
    % dg(B)/dRf
    % this matches the numerical solution, but I don't know why.
    dgB_dRf = zeros(N*M, 1);
    for i = 1:N
        for j = 1:M
            dgB_dRf(N*(j-1)+i, 1) = sum(-(1 - Beta(i,:)*iotaMp)*(Rmp(j,:) - iotaMp(j)*zbrate) ...
                 - e(i, :)*iotaMp(j));
        end
    end
    dgB_dRf = dgB_dRf/T;

    % dg(B)/dgamma
    dgB_dgamma = zeros(N*M, K);

    for i = 1:N
        for j = 1:M
            for k = 1:K
                dgB_dgamma(N*(j-1)+i, k) = sum(-(1 - Beta(i,:)*iotaMp)*Z(k, :)*(Rmp(j,:) - iotaMp(j)*zbrate)') ...
                     - e(i, :)*iotaMp(j)*Z(k, :)';
            end
        end
    end

    dgB_dgamma = dgB_dgamma/T;


    % dg(B)/drho
    dgB_drho = zeros(N*M,1);



    % dg(B)/dB
    dgB_dB = zeros(N*M, N*M);
    
    for i = 1:N
        for j = 1:M
            for p = 1:N
                for q = 1:M
                    % dgB_dB(N*(j-1)+i, N*(q-1)+p) = Sinv(i,p)*(iotaM(q)*zbrate-Rm(q,:))*(Rm(j,:)-iotaM(j)*zbrate)';
                    % No need to include the covariance matrix?
                    if i == p
                        dgB_dB(N*(j-1)+i, N*(q-1)+p) = (iotaMp(q)*zbrate-Rmp(q,:))*(Rmp(j,:)-iotaMp(j)*zbrate)';
                    end
                end
            end
        end
    end
    
    dgB_dB = dgB_dB/T;

    % dg(B)/dalpha
    dgB_dalpha = zeros(N*M, N);

    for i = 1:N
        for j = 1:M
            for p = 1:N
                if i == p
                    dgB_dalpha(N*(j-1)+i, p) = sum(-(Rmp(j,:) - iotaMp(j)*zbrate));
                end
            end
        end
    end
    
    dgB_dalpha = dgB_dalpha/T;

    momjacobian = [dge_dRf, dge_dgamma, dge_drho, dge_dB, dge_dalpha;...
                   dgB_dRf, dgB_dgamma, dgB_drho, dgB_dB, dgB_dalpha];

    