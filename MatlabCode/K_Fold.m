function [errors] = K_Fold(Rinput, Rminput, Zinput, Rbinput, iotaN, iotaM, ConsG,inflation,Rfex,psi,sig,asset,har,NLConsFactor,SigmaType,Folds)

N = size(Rinput, 1); T = size(Rinput, 2); M = size(Rminput, 1); K = size(Zinput, 1);

errors = zeros(Folds,1);


batch_size = floor(T/Folds);
indx = 1:T;
    
for i = 1:Folds
    if i ~= Folds
        test_idx = indx(batch_size*(i-1)+1: batch_size*i);
    else 
        test_idx = indx(batch_size*(i-1)+1: end);
    end
    train_idx = indx(~ismember(indx, test_idx));

    [~,Theta1,~,~,weight] = InstrumentGMMWrapperConc(Rinput(:,train_idx), Rminput(:,train_idx), Zinput(:,train_idx), Rbinput(:, train_idx), iotaN, iotaM, ConsG(train_idx),inflation(train_idx),Rfex(train_idx),psi,sig,asset,har,NLConsFactor,SigmaType);


    Rf = Theta1(1);
    gamma = reshape(Theta1(2:end-2), K, 1);

    zbrate = gamma'*Zinput(:,test_idx) + Rf +  Rbinput(:,test_idx);
    portRet = weight *  Rinput(:,test_idx);

    errors(i) = 1/batch_size*(portRet-zbrate)*(portRet-zbrate)';

end



end



