%test code to validate analytical portions of Jacobians
%to run this code, you must add in 'jac' as a third argument to return in
%InstGMMSE. it is not enabled by default for performance reasons.

load("../Output/MainRun_Main.mat");

%test jacobian for zero-beta rate
[pvar,mvar, jac] = InstGMMSE(Theta,Beta,alphas, Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, cons_gr_ann/12,inflation,weight,opts.har,opts.SigmaInit,opts.NLConsFactor);


 xinit = [Theta(1:end-1);reshape(Beta,size(Beta,1)*size(Beta,2),1);alphas'];

 func = @(x) InstMomentsConc([x(1:2+K);opts.SigmaInit], ... %extract thetas
     reshape(x(3+K+size(Beta,1)*size(Beta,2):end),size(Beta,1),1)', ... %extract alphas
     reshape(x(3+K:2+K+size(Beta,1)*size(Beta,2)),size(Beta,1),size(Beta,2)), ... %extract betas
     weight, Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, cons_gr_ann/12,inflation,[],'ZBFull',opts.NLConsFactor);
 
jac2 = jacobianest(func,xinit);

diffj = abs(jac-jac2);
disp('max abs. error between jacobians')
max(diffj,[],"all")


%test jacobian for risk-free rate 

[~,ThetaRF,sviRF,~, ~, ~, ~, ~, mmomentsRF] = InstrumentGMMWrapperConc(Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, cons_gr_ann/12, inflation,Rfex,0,opts.SigmaInit,'RF',opts.har,opts.NLConsFactor);

 funcRF = @(x) InstMomentsConc([x;opts.SigmaInit], ... %extract thetas
     [], ... %extract alphas
     [], ... %extract betas
     weight, Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, cons_gr_ann/12,inflation,Rfex,'RF',opts.NLConsFactor);
 
jac2RF = jacobianest(funcRF,ThetaRF(1));

%code from analytical computation
Nm = 1;
e1 = [1;zeros(K,1)];
jacRF = -100*(e1+mmomentsRF(Nm:end)/100)';

diffjRF = abs(jacRF'-jac2RF);
disp('max abs. error between RF jacobians')
max(diffjRF,[],"all")