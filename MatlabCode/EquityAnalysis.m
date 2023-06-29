load("../Output/MainRun_AltDP.mat");


rho = (0.96)^(1/12);
sigma = 5;

% univariate analysis
lm_zb = fitlm(log(1+zbrateReal(1:end-1)/100),log(1+zbrateReal(2:end)/100));

phi_zb = lm_zb.Coefficients.Estimate(2);
sigma_zb = lm_zb.RMSE;

%annualized, in percent
univol_zb = sigma_zb / (1-rho*phi_zb) * 100 * sqrt(12) * (1-1/sigma)


% same for treasury
lm_rf = fitlm(log(1+RbinputReal(1:end-1)/100),log(1+RbinputReal(2:end)/100));

phi_rf = lm_rf.Coefficients.Estimate(2);
sigma_rf = lm_rf.RMSE;

%annualized, in percent
univol_rf = sigma_rf / (1-rho*phi_rf) * 100 * sqrt(12) * (1-1/sigma)




%VAR analysis

dpr = Instruments{:,'DP_ratio'};
mrets = table2array(Factors(:,'Mkt'));


divg = 100*log( 12*(1+mrets(2:end)/100)./dpr(1:end-1)./(1+1./dpr(2:end)*12)) - cpi(2:end);

divg_cs = 100*( log(1+mrets(2:end)/100) + log(rho) + (1-rho)*log(1/rho-1) + rho * log(dpr(2:end)/12) - log(dpr(1:end-1)/12)) - cpi(2:end);

mbill = movmean(Rbinput,[11 0]);
rrel = 12*(Rbinput(13:end) - mbill(12:end-1));

%original campbell specification with real returns
Z0 = [Rminput(1,12:end-1)-cpi(2+12:end-1)';Zinput(4,13:end);rrel];

erm0 = [1;0;0];
Phi0=(Z0(:,1:end-1)' \ Z0(:,2:end)');
resid0 = Z0(:,2:end)' - Z0(:,1:end-1)'*Phi0;
Sigma_mat0 = cov(resid0);
imat0 = inv(eye(length(erm0))-rho*Phi0');

lam_rm0 = rho*erm0' * Phi0'*imat0;
lam_drf0 = erm0' + lam_rm0;
vec_0 = [lam_rm0; lam_drf0];
var_rf0 = vec_0 * Sigma_mat0 * vec_0' / Sigma_mat0(1,1);

results_0 = [1-Sigma_mat0(1,1)/var(Z0(1,:)), var_rf0(2,2), var_rf0(1,1),-2*var_rf0(1,2)]


%original campbell specification with excess returns

Z1 = [RbinputReal(13:end);Zinput(4,13:end);rrel; Rminput(1,12:end-1)-Rbinput(12:end-1)];

ef1 = [1;0;0;0];
erm1 = [0;0;0;1];


Phi1=(Z1(:,1:end-1)' \ Z1(:,2:end)');
resid1 = Z1(:,2:end)' - Z1(:,1:end-1)'*Phi1;
Sigma_mat1 = cov(resid1);
imat1 = inv(eye(length(ef1))-rho*Phi1');

lam_rf1 = ef1' * imat1;
lam_rm_rf1 = rho*erm1' * Phi1'*imat1;
lam_drf1 = erm1' + lam_rf1 + lam_rm_rf1;
vec_rf1 = [lam_rf1; lam_rm_rf1; lam_drf1];

var_rf1 = vec_rf1 * Sigma_mat1 * vec_rf1' / Sigma_mat1(4,4);

results_1 = [1-Sigma_mat1(4,4)/var(Z1(4,:)), var_rf1(3,3), var_rf1(1,1), var_rf1(2,2), -2*var_rf1(3,1), -2*var_rf1(3,2),2*var_rf1(1,2)]

sum(results_1(2:end))


Z2 = [Zinput(:,2:end);Rminput(1,1:end-1)-Rbinput(1:end-1)];

inf_vec = [100*beta_inf(1:end-1).*Zscales';0];

gamma2 = [gamma;0];
ef = zeros(size(gamma2));
ef(1) = Zscales(1);
ef = ef + inf_vec;

erm_rf = zeros(size(gamma2));
erm_rf(end) = 1;

erm_zb = erm_rf - gamma2;


Phi=(Z2(:,1:end-1)' \ Z2(:,2:end)');

resid = Z2(:,2:end)' - Z2(:,1:end-1)'*Phi;

Sigma_mat = cov(resid);


imat = inv(eye(length(gamma2))-rho*Phi');

lam_zb = (ef+gamma2)' * imat;

lam_rf = ef' * imat;

lam_rm_zb = rho*erm_zb' * Phi'*imat;
lam_rm_rf = rho*erm_rf' * Phi'*imat;

lam_dzb = erm_zb' + lam_zb + lam_rm_zb;
lam_drf = erm_rf' + lam_rf + lam_rm_rf;

vec_zb = [lam_zb; lam_rm_zb; lam_dzb];
vec_rf = [lam_rf; lam_rm_rf; lam_drf];

var_zb = vec_zb * Sigma_mat * vec_zb' / Sigma_mat(end,end);

var_rf = vec_rf * Sigma_mat * vec_rf' / Sigma_mat(end,end);

results_zb = [1-Sigma_mat(end,end)/var(Z2(end,:)), var_zb(3,3), var_zb(1,1), var_zb(2,2), -2*var_zb(3,1), -2*var_zb(3,2),2*var_zb(1,2)]

results_rf = [1-Sigma_mat(end,end)/var(Z2(end,:)), var_rf(3,3), var_rf(1,1), var_rf(2,2), -2*var_rf(3,1), -2*var_rf(3,2),2*var_rf(1,2)]

% varvol_zb = sqrt((ef+gamma2)' * imat*Sigma_mat*imat' * (ef+gamma2))* sqrt(12)
% 
% varvol_rf = sqrt(ef' * imat*Sigma_mat*imat' * ef)* sqrt(12) 
% 
% varvol_erm = sqrt(erm' * imat*Sigma_mat*imat' * erm)* sqrt(12) 
% 
% 
% dvol_zb = sqrt((erm+((erm+ef+gamma2)' * imat)*Sigma_mat*imat' * (ef+gamma2))* sqrt(12)