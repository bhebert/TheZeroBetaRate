function [errors,errors_c] = K_Fold(Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, ConsG,inflation,Rfex,psi,psic,sig,asset,har,NLConsFactor,Folds,RandomFeature)

N = size(Rinput, 1); T = size(Rinput, 2); M = size(Rminput, 1); K = size(Zinput, 1);

errors = zeros(Folds,1);
errors_c = zeros(Folds,1);

batch_size = floor(T/Folds);
indx = 1:T;
    
for i = 1:Folds
    if i ~= Folds
        test_idx = indx(batch_size*(i-1)+1: batch_size*i);
    else 
        test_idx = indx(batch_size*(i-1)+1: end);
    end
    train_idx = indx(~ismember(indx, test_idx));

    if RandomFeature
        [Theta1,~,~,~,weight,cbetas] = RandomFeatureEst(Rinput(:,train_idx), Rminput(:,train_idx), Zinput(:,train_idx), Rbinput(:, train_idx),ConsG(train_idx),inflation(train_idx),iotaN, iotaM, sig,NLConsFactor,psi,psic);
        p_cons = cbetas(1) + cbetas(2:end)'*Zinput(:,test_idx);
    else
        [~,Theta1,~,~,weight] = InstrumentGMMWrapperConc(Rinput(:,train_idx), Rminput(:,train_idx), Zinput(:,train_idx), Rbinput(:, train_idx), iotaN, iotaM, ConsG(train_idx),inflation(train_idx),Rfex(train_idx),psi,sig,asset,har,NLConsFactor);
        p_cons = ConsG(test_idx)';
    end

    Rf = Theta1(1);
    gamma = reshape(Theta1(2:end-2), K, 1);

    zbrate = gamma'*Zinput(:,test_idx) + Rf +  Rbinput(:,test_idx);
    portRet = weight *  Rinput(:,test_idx);

     

    errors(i) = (1/batch_size*(portRet-zbrate)*(portRet-zbrate)') / var(portRet);
    errors_c(i) = (1/batch_size*(ConsG(test_idx)'-p_cons)*(ConsG(test_idx)-p_cons')) / var(ConsG(test_idx));
end



end



