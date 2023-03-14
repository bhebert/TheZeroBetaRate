clear all
close all
clc

%% Parameters to play with

global sigma nu psi

sigma = 5;
nu = 1/6;
psi = 10;
M(1) = 1;
M(2) = 0.99;
B(1) =1.02;
B_ss_annual = .9;
M_ss_annual = .3;


%% Fixed Parameters
delta_annual = 0.92;
mubar_annual = 0.02;
rho_annual = 0.1;
spread_annual = 0.07;
% sticky_annual = 0.8;
% sticky_bond_annual = 0.5;
T_annual = 300;
Tplot_annual = 4;
c_ss_annual = 1;


freq = 1;

delta = delta_annual^(1/freq);
mubar = (1+mubar_annual)^(1/freq)-1;
rho = rho_annual^(1/freq);
spread = (1+spread_annual)^(1/freq)-1;
% sticky = sticky_annual^(1/freq);
% sticky_bond = sticky_bond_annual^(1/freq);
T = T_annual*freq;
B_ss = B_ss_annual;
M_ss = M_ss_annual;
Tplot = Tplot_annual*freq;
c_ss = c_ss_annual/freq;
ts = linspace(1,Tplot,Tplot);

P_ss = 1;
r_ss = 1/delta*(1+mubar)-1;
rb_ss = r_ss - spread*(1+r_ss);

guesspar = [.25;.1];

[par err] = fsolve(@(e) errss_alt(e,P_ss,c_ss,M_ss,B_ss,r_ss,rb_ss),guesspar)
                               
alpha = par(1);
eta = par(2);


%% Path of shocks

P(1) = P_ss;
B(1) = B(1)*B_ss;

for j = 3:T
    M(j) = M(j-1);
end

M = M*M_ss.*(1+mubar).^(0:T-1);

H_ss = (alpha^(1/nu) * M_ss^((nu-1)/nu) + (1-alpha)^(1/nu) * B_ss^((nu-1)/nu))^(nu/(nu-1));

iota_ss = r_ss/(1+r_ss)*alpha^(-1/nu)*(M_ss/H_ss)^(1/nu);

H(1) = (alpha^(1/nu) * M(1)^((nu-1)/nu) + (1-alpha)^(1/nu) * B(1)^((nu-1)/nu))^(nu/(nu-1));

%% Solve forward

dist = 10;
rL = 0.0001;
rH = 0.2;
explode = 0.2/freq;
e = 1e-5;
k = 1;
maxiter = 150;

while (dist > e)*(k<maxiter)
    r(1) = (rL +rH)/2;
    iota(1) = r(1)/(1+r(1))*alpha^(-1/nu)*(M(1)/H(1))^(1/nu);
    s(1) = iota(1)*(B(1)/H(1))^(-1/nu)*(1-alpha)^(1/nu);
    c(1) = eta^(-1/sigma)*(H(1)/P(1))^(psi/sigma)*(iota(1))^(1/sigma);
    rb(1) = r(1) - s(1)*(1+r(1));
    
    for t = 2:T
        B(t) = M(t)*B_ss/M_ss;
        P(t) = P_ss*(1+mubar)^(t-1);
        H(t) = (alpha^(1/nu) * M(t)^((nu-1)/nu) + (1-alpha)^(1/nu) * B(t)^((nu-1)/nu))^(nu/(nu-1));
        r(t) = delta*r(t-1)/...
            ( (H(t)/H(t-1))^(psi-1/nu)*(M(t)/M(t-1))^(1/nu)*(P(t)/P(t-1))^(1-psi) - delta*r(t-1));	
        iota(t) = r(t)/(1+r(t))*alpha^(-1/nu)*(M(t)/H(t))^(1/nu);
        s(t) = iota(t)*(B(t)/H(t))^(-1/nu)*(1-alpha)^(1/nu);
        c(t) = eta^(-1/sigma)*(H(t)/P(t))^(psi/sigma)*(iota(t))^(1/sigma);
        rb(t) = r(t) - s(t)*(1+r(t));


    end
    
    if real(r(T))>r_ss
        rH = r(1);
    elseif max(real(r)) > explode
        rH = r(1);
    else
        rL = r(1);
    end
       
    if max(real(r)) > explode
        dist = 10;
    else
        dist = (real(r(T)) - r_ss)^2;
    end
    
     k = k+1;
end  

% Compute inflation and real interest rates
pi = P(1:T)./[1/(1+mubar),P(1:T-1)]-1;
r_real = r-pi;

 
%% Plot 

figure(2)

% subplot(2,5,1)
% plot(ts,(1+mu(1:Tplot)).^freq-1,ts,(1+mubar)^freq*ones(1,Tplot)-1,'--')
% xlabel('t')
% ylabel('\mu')

subplot(2,4,2)
plot((ts-1)/freq,(1+r(1:Tplot)).^freq-1,(ts-1)/freq,(1+r_ss)^freq*ones(1,Tplot)-1,'--')
xlabel('t')
ylabel('r')

subplot(2,4,1)
plot((ts-1)/freq,M(1:Tplot)./(1+mubar).^(0:Tplot-1)/M_ss)
xlabel('t')
ylabel('M')

% subplot(2,5,5)
% plot(ts,P(1:Tplot)./(1+mubar).^(0:Tplot-1))
% xlabel('t')
% ylabel('P')
% 
subplot(2,4,4)
plot((ts-1)/freq,c(1:Tplot)*freq)
xlabel('t')
ylabel('c')
% % 
% subplot(2,5,6)
% plot(ts,(1+pi(1:Tplot)).^freq-1)
% xlabel('t')
% ylabel('\pi')
% 
subplot(2,4,6)
plot((ts-1)/freq,(1+s(1:Tplot)).^freq-1, (ts-1)/freq,spread_annual*ones(1,Tplot),'--')
xlabel('t')
ylabel('spread')
% 
subplot(2,4,3)
plot((ts-1)/freq,(1+rb(1:Tplot)).^freq-1,(ts-1)/freq,(1+rb_ss)^freq*ones(1,Tplot)-1,'--')
xlabel('t')
ylabel('r_b')


subplot(2,4,5)
plot((ts-1)/freq,B(1:Tplot)./(1+mubar).^(0:Tplot-1)/B_ss)
xlabel('t')
ylabel('Bond supply')

subplot(2,4,7)
plot((ts-1)/freq,H(1:Tplot)./(1+mubar).^(0:Tplot-1)/H_ss)
xlabel('t')
ylabel('Composite Liquidity Supply')


subplot(2,4,8)
plot((ts-1)/freq,(1+iota(1:Tplot)).^freq-1)
xlabel('t')
ylabel('Composite Liquidity Cost')

saveas(gcf, "../Output/model.png")

%% On-impact:

% This is a RR shock
% dr0 = [((1+r(1))^freq - (1+r_ss)^freq)*100 ,-1.514] 
% drb = [((1+rb(1))^freq - (1+rb_ss)^freq)*100, 0.338]
% dgc = [((c(2)/c(1))^freq - 1)*100, -.292]

% This is 1/10 of a NS shock
dr0 = [((1+r(1))^freq - (1+r_ss)^freq)*100 ,-2.5*5/10] 
drb = [((1+rb(1))^freq - (1+rb_ss)^freq)*100, 1.8/10]
dgc = [((c(2)/c(1))^freq - 1)*100, -.2]


% dlogc = [(log(c(1))-log(c_ss))*100, -5.95] 

return
