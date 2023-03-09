clear;
clc;

Nlags = 25;

%% Load data
load("../Output/MainRun_Ridge.mat");
Zadj = gamma .* Zinput;
p_cons = p_cons';
ZBR = table(dts', zbrateReal', RbinputReal',Zadj',p_cons');
% ZBR = table(dts, riskfree', Rf(1:end-1)*12);
ZBR = renamevars(ZBR, ["Var1", "Var2", "Var3","Var5"], ["Date", "ZBR", "RF","PCons"]);

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
for l = 0:Nlags-1
    % ---- zero beta rate ----
    % change
    %
    mdl = fitlm(REGDT(2:end-l-1,6),12*(REGDT(l+3:end,1) - REGDT(1:end-2-l,1)));

    reg_results(1+l, 1, 1) = mdl.Coefficients.Estimate(2);
    CIs = coefCI(mdl);
    reg_results(1+l, 2, 1) = CIs(2,1); 
    reg_results(1+l, 3, 1) = CIs(2,2); 

    % ---- t-bill rate ----
    mdl = fitlm(REGDT(2:end-l-1,6),12*(REGDT(l+3:end,2) - REGDT(1:end-2-l,2)));

    reg_results(1+l, 1, 2) = mdl.Coefficients.Estimate(2);
    CIs = coefCI(mdl);
    reg_results(1+l, 2, 2) = CIs(2,1); 
    reg_results(1+l, 3, 2) = CIs(2,2); 

    % ---- exp. cons ----
    mdl = fitlm(REGDT(2:end-l-1,6),(REGDT(l+3:end,3) - REGDT(1:end-2-l,3)));

    reg_results(1+l, 1, 3) = mdl.Coefficients.Estimate(2);
    CIs = coefCI(mdl);
    reg_results(1+l, 2, 3) = CIs(2,1); 
    reg_results(1+l, 3, 3) = CIs(2,2); 
end

figure(1)
% zero beta rate
tiledlayout(1,2)
ax1 = nexttile;
hold on
plot(0:Nlags-1, reg_results(:, 1, 1), "k-")
plot(0:Nlags-1, reg_results(:, 2, 1), "b--")
plot(0:Nlags-1, reg_results(:, 3, 1), "b--")
plot(0:Nlags-1, reg_results(:, 1, 2), "r-")
plot(0:Nlags-1, reg_results(:, 2, 2), "g--")
plot(0:Nlags-1, reg_results(:, 3, 2), "g--")
hold off
yline(0, 'r-')
xlabel('Lag')
ylabel('Coefficient')
title("Zero Beta Rate and T-Bill Yield")
pbaspect(ax1,[1 1 1])

% t-bill rate
ax2 = nexttile;
hold on
plot(0:Nlags-1, reg_results(:, 1, 3), "k-")
plot(0:Nlags-1, reg_results(:, 2, 3), "b--")
plot(0:Nlags-1, reg_results(:, 3, 3), "b--")
hold off
yline(0, 'r-')
xlabel('Lag')
ylabel('Coefficient')
title("Convenience Yield")
pbaspect(ax2,[1 1 1])




saveas(gcf, "../Output/local_projection_rr.png")

% Local Projections With NS
REGDT = table2array(DT2(:, 2:end));
reg_results = zeros(Nlags, 3, 2);
for l = 0:Nlags-1
    % ---- zero beta rate ----
    % change
    mdl = fitlm(REGDT(2:end-l-1,5),12*(REGDT(l+3:end,1) - REGDT(1:end-2-l,1)));

    reg_results(1+l, 1, 1) = mdl.Coefficients.Estimate(2);
    CIs = coefCI(mdl);
    reg_results(1+l, 2, 1) = CIs(2,1); 
    reg_results(1+l, 3, 1) = CIs(2,2); 

    % ---- t-bill rate ----


    mdl = fitlm(REGDT(2:end-l-1,5),12*(REGDT(l+3:end,2) - REGDT(1:end-2-l,2)));


    reg_results(1+l, 1, 2) = mdl.Coefficients.Estimate(2);
    CIs = coefCI(mdl);
    reg_results(1+l, 2, 2) = CIs(2,1); 
    reg_results(1+l, 3, 2) = CIs(2,2); 


        % ---- exp. cons ----
    mdl = fitlm(REGDT(2:end-l-1,5),(REGDT(l+3:end,3) - REGDT(1:end-2-l,3)));

    reg_results(1+l, 1, 3) = mdl.Coefficients.Estimate(2);
    CIs = coefCI(mdl);
    reg_results(1+l, 2, 3) = CIs(2,1); 
    reg_results(1+l, 3, 3) = CIs(2,2); 
end

figure(2)


% zero beta rate
tiledlayout(1,2)
ax1 = nexttile;
hold on
plot(0:Nlags-1, reg_results(:, 1, 1), "k-")
plot(0:Nlags-1, reg_results(:, 2, 1), "b--")
plot(0:Nlags-1, reg_results(:, 3, 1), "b--")
plot(0:Nlags-1, reg_results(:, 1, 2), "r-")
plot(0:Nlags-1, reg_results(:, 2, 2), "g--")
plot(0:Nlags-1, reg_results(:, 3, 2), "g--")
hold off
yline(0, 'r-')
xlabel('Lag')
ylabel('Coefficient')
title("Zero Beta Rate and T-Bill Yield")
pbaspect(ax1,[1 1 1])

% t-bill rate
ax2 = nexttile;
hold on
plot(0:Nlags-1, reg_results(:, 1, 3), "k-")
plot(0:Nlags-1, reg_results(:, 2, 3), "b--")
plot(0:Nlags-1, reg_results(:, 3, 3), "b--")
hold off
yline(0, 'r-')
xlabel('Lag')
ylabel('Coefficient')
title("Convenience Yield")
pbaspect(ax2,[1 1 1])


% % zero beta rate
% tiledlayout(1,2)
% ax1 = nexttile;
% hold on
% plot(0:Nlags-1, reg_results(:, 1, 1), "k-")
% plot(0:Nlags-1, reg_results(:, 2, 1), "b--")
% plot(0:Nlags-1, reg_results(:, 3, 1), "b--")
% hold off
% yline(0, 'r-')
% xlabel('Lag')
% ylabel('Coefficient')
% title("Zero Beta Rate")
% pbaspect(ax1,[1 1 1])
% 
% % t-bill rate
% ax2 = nexttile;
% hold on
% plot(0:Nlags-1, reg_results(:, 1, 2), "k-")
% plot(0:Nlags-1, reg_results(:, 2, 2), "b--")
% plot(0:Nlags-1, reg_results(:, 3, 2), "b--")
% hold off
% yline(0, 'r-')
% xlabel('Lag')
% ylabel('Coefficient')
% title("T-Bill Rate")
% pbaspect(ax2,[1 1 1])

saveas(gcf, "../Output/local_projection_ns.png")


