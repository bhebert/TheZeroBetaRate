clear;
clc;
close all;

% color list
% visually distinct colors via colorbrewer2.org
colors_list = {'#e41a1c','#377eb8', '#4daf4a','#984ea3','#ff7f00'};

%flag for including contemporaneous result
%0: don't use contemporaneous (WP version)
%1: include comptemporaneous (suggested by referee)
use_0 = 1;

%flag for using specification with lagCPI
use_lagcpi = false; %false: main text true: appendix


Nlags = 12;

%% Load data

if use_lagcpi
    load("../Output/MainRun_AltCPI.mat");
    ext='_altcpi';
else
    load("../Output/MainRun_Main.mat");
    ext='';
end
%


%fix due to redefinition of gamma in paper
Zadj = (gammaNC + [1/ZSDiagNC(2,2); zeros(K-1,1)]) .* ZinputNC;


ZBR = table(dtsNC', zbrateRealNC', RbinputRealNC',Zadj',p_consNC',zbrateNC');
% ZBR = table(dts, riskfree', Rf(1:end-1)*12);
ZBR = renamevars(ZBR, ["Var1", "Var2", "Var3","Var5","Var6"], ["Date", "ZBR", "RF","PCons","ZBN"]);

RR = readtable("../Input/RR_shocks.csv");
NS = readtable("../Input/NS_shocks.csv");

ZBR2 = ZBR(:,["Date","ZBR","RF","PCons"]);

DT1 = innerjoin(ZBR2, RR, "Keys",'Date');
DT2 = innerjoin(ZBR2, NS, "Keys",'Date');


DT1a = innerjoin(ZBR, RR, "Keys",'Date');
DT2a = innerjoin(ZBR, NS, "Keys",'Date');

writetable(DT1a,'../Output/RRData_Main.csv');
writetable(DT2a,'../Output/NSData_Main.csv');

%% Local Projections With RR
REGDT = table2array(DT1(:, 2:end));
reg_results = zeros(Nlags, 3, 3);
for l = 0:Nlags-1+use_0
    % ---- zero beta rate ----
    % change
    %
    mdl = fitlm(REGDT(2:end-l-1,6),12*(REGDT(l+3-use_0:end-use_0,1) - REGDT(1:end-2-l,1)));

    reg_results(1+l, 1, 1) = mdl.Coefficients.Estimate(2);
    CIs = coefCI(mdl);
    reg_results(1+l, 2, 1) = CIs(2,1); 
    reg_results(1+l, 3, 1) = CIs(2,2); 

    % ---- t-bill rate ----
    mdl = fitlm(REGDT(2:end-l-1,6),12*(REGDT(l+3-use_0:end-use_0,2) - REGDT(1:end-2-l,2)));

    reg_results(1+l, 1, 2) = mdl.Coefficients.Estimate(2);
    CIs = coefCI(mdl);
    reg_results(1+l, 2, 2) = CIs(2,1); 
    reg_results(1+l, 3, 2) = CIs(2,2); 

    % ---- exp. cons ----
    mdl = fitlm(REGDT(2:end-l-1,6),(REGDT(l+3-use_0:end-use_0,3) - REGDT(1:end-2-l,3)));

    reg_results(1+l, 1, 3) = mdl.Coefficients.Estimate(2);
    CIs = coefCI(mdl);
    reg_results(1+l, 2, 3) = CIs(2,1); 
    reg_results(1+l, 3, 3) = CIs(2,2); 
end

testSigma = 10;

figure(1);
% zero beta rate

hold on

plot((1:Nlags+use_0)-use_0, (reg_results(:, 1, 1))/testSigma, '-','LineWidth',2,'Color',colors_list{1})
plot((1:Nlags+use_0)-use_0, (reg_results(:, 1, 2)), '-.','LineWidth',2,'Color',colors_list{3})
plot((1:Nlags+use_0)-use_0, (reg_results(:, 1, 3)), '--','LineWidth',2,'Color',colors_list{2})

set(gca,'TickLabelInterpreter','latex')
legend('Real Zero-Beta Rate /'+sprintf("%d",testSigma),'$\bf{E}[$Real T-Bill Return$]$','$\bf{E}[$Cons. Growth$]$','show','Location','southwest','Interpreter','Latex')

hold off
yline(0, 'k-','HandleVisibility','off');
xlabel('Months after Shock','Interpreter','latex')
ylabel('Percent','Interpreter','latex')
title("Effects of the Romer-Romer Shock",'Interpreter','latex')
ylim([-6,4])
saveas(gcf, "../Output/local_projection_rr"+ext+".png")

% Local Projections With NS
REGDT = table2array(DT2(:, 2:end));
reg_results = zeros(Nlags, 3, 2);
for l = 0:Nlags-1+use_0
    % ---- zero beta rate ----
    % change
    mdl = fitlm(REGDT(2:end-l-1,5),12*(REGDT(l+3-use_0:end-use_0,1) - REGDT(1:end-2-l,1)));

    reg_results(1+l, 1, 1) = mdl.Coefficients.Estimate(2);
    CIs = coefCI(mdl);
    reg_results(1+l, 2, 1) = CIs(2,1); 
    reg_results(1+l, 3, 1) = CIs(2,2); 

    % ---- t-bill rate ----


    mdl = fitlm(REGDT(2:end-l-1,5),12*(REGDT(l+3-use_0:end-use_0,2) - REGDT(1:end-2-l,2)));


    reg_results(1+l, 1, 2) = mdl.Coefficients.Estimate(2);
    CIs = coefCI(mdl);
    reg_results(1+l, 2, 2) = CIs(2,1); 
    reg_results(1+l, 3, 2) = CIs(2,2); 


        % ---- exp. cons ----
    mdl = fitlm(REGDT(2:end-l-1,5),(REGDT(l+3-use_0:end-use_0,3) - REGDT(1:end-2-l,3)));

    reg_results(1+l, 1, 3) = mdl.Coefficients.Estimate(2);
    CIs = coefCI(mdl);
    reg_results(1+l, 2, 3) = CIs(2,1); 
    reg_results(1+l, 3, 3) = CIs(2,2); 
end

figure(2)


hold on
plot((1:Nlags+use_0)-use_0, (reg_results(:, 1, 1))/testSigma, '-', 'LineWidth',2,'Color',colors_list{1})
plot((1:Nlags+use_0)-use_0, (reg_results(:, 1, 2)), '-.', 'LineWidth',2,'Color',colors_list{3})
plot((1:Nlags+use_0)-use_0, (reg_results(:, 1, 3)), '--', 'LineWidth',2,'Color',colors_list{2})
ylim([-6,4])
set(gca,'TickLabelInterpreter','latex')
legend('Real Zero-Beta Rate /'+sprintf("%d",testSigma),'$\bf{E}[$Real T-Bill Return$]$','$\bf{E}[$Cons. Growth$]$','show','Location','southwest','Interpreter','Latex')

yline(0, 'k-','HandleVisibility','off');
hold off
xlabel('Months after Shock','Interpreter','latex')
ylabel('Percent','Interpreter','latex')
title("Effects of the Nakamura-Steinsson Shock",'Interpreter','latex')


saveas(gcf, "../Output/local_projection_ns"+ext+".png")


