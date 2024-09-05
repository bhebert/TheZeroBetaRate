clear
close all
clc

global sigma nu psi

colors_list = {'#e41a1c','#377eb8', '#4daf4a','#984ea3','#ff7f00'};

%% Parameters

sigma = 7.5;

delta = 1/1.1;                  % Discount factor
mubar = 0.02;                   % Steady state inflation
spread = 0.08;                  % Steady state spread

%% Precision parameters
dist = 10;          
e = 1e-10;                 % Convergence criterion
T = 150;                   % Horizon of simulation
low = 0.9;                 % Bounds of initial guess
high = 1.01;
explode = 10;               % Detect an explosive path

%% Targets for steady state
B_ss = .5/.68;                  % Target: public debt to consumption ratio
M_ss = .16/.68;                 % Target: M1 to consumption ratio
c_ss = 1;                       % Normalization
P_ss = 1;                       % Normalization
r_ss = 1/delta*(1+mubar)-1;
rb_ss = r_ss - spread*(1+r_ss);

%% Finde etaM and etaB that match steady state targets

etaM = M_ss*c_ss^-sigma*r_ss;
etaB = B_ss*c_ss^-sigma*spread;

%% Path of shocks

% Paper revision with sigma = 7.5 version
M(1) = 1;           
M(2) = .9819;      
B(1) = 1.0073;     

M(3:T) = M(2);      % M stays constant after period 2

P(1) = P_ss;
B(1) = B(1)*B_ss;   % Rescale by steady state bond supply

M = M*M_ss.*(1+mubar).^(0:T-1);     % M scales up by steady state inflation

%% Solve forward. Guess a value for the zero-beta rate in period 1, solve forward and update

cL = c_ss*low;
cH = c_ss*high;

while dist > e

    % Update guess of initial consumption
    c(1) = (cL + cH)/2;

    % Zero-beta rate from money demand
    r(1) = etaM/(M(1)/P(1))*c(1)^sigma;

    % Compute spread from bond demand
    s(1) = etaB/(B(1)/P(1))*c(1)^sigma;

    % Compute safe rate from r0 and the spread
    rb(1) = r(1) - s(1)*(1+r(1));
        
    for t = 2:T

        % Price level fixed (follows steady state inflation)
        P(t) = P_ss*(1+mubar)^(t-1);

        % Consumption from Euler equation:
        c(t) = c(t-1)*(delta*(1+r(t-1))/(1+mubar))^(1/sigma);

        % Zero-beta rate from money demand
        r(t) = etaM/(M(t)/P(t))*c(t)^sigma;

        % Spread goes to steady state
        s(t) = spread;

        % Bond supply to match spread
        B(t) = etaB*c(t)^sigma/s(t)*(1+mubar)^(t-1);

        % Safe rate from spread 
        rb(t) = r(t) - s(t)*(1+r(t));

       
    end

    
    if real(c(T)) > (M(2)/(M_ss*(1+mubar)))^(1/sigma)
        cH = c(1);
    elseif max(real(c)) > explode
        cH = c(1);
    else
        cL = c(1);
    end
     
    % Compute distance to convergence
    if max(real(r)) > explode        
        dist = explode;            
    else
        dist = (real(r(T)) - r_ss)^2;
    end
   
   
end  

%% Copmute inflation and real rates

pi = [mubar,P(2:T)./P(1:T-1)-1];
r_real = (1+r)./(1+pi)-1;
rb_real = (1+rb)./(1+pi)-1;

 
%% Compute on-impact effect of a shock

% Zero-beta rate
dr0 = ((1+r(1)) - (1+r_ss))*100

% Safe rate
drb = ((1+rb(1)) - (1+rb_ss))*100

% Consumption level
dlogc = (log(c(1))-log(c_ss))*100 

% Consumption growth
dgc = (c(2)/c(1)-1)*100 


%% 3-period plots

Tplot = 3;

ts = linspace(0,Tplot,Tplot+1);

f = figure(1)

subplot(2,3,1)
plot(ts,[M_ss/(1+mubar),M(1:Tplot)]./(1+mubar).^(ts-1)/M_ss,'Color',colors_list{2},'LineWidth',2)
xlabel('Year','Interpreter','Latex')
ylabel('Money Supply','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')


subplot(2,3,4)
plot(ts,[B_ss/(1+mubar),B(1:Tplot)]./(1+mubar).^(ts-1)/B_ss,'Color',colors_list{2},'LineWidth',2)
xlabel('Year','Interpreter','Latex')
ylabel('Bond Supply','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')


subplot(2,3,3)
plot(ts,(1+[(1+r_ss)/(1+mubar)-1,r_real(1:Tplot)])-1,'Color',colors_list{2},'LineWidth',2)
xlabel('Year','Interpreter','Latex')
ylabel('Real Zero-Beta Rate','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')


subplot(2,3,6)
plot(ts,(1+[(1+rb_ss)/(1+mubar)-1,rb_real(1:Tplot)])-1,'Color',colors_list{2},'LineWidth',2)
xlabel('Year','Interpreter','Latex')
ylabel('Real Safe Rate','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')


subplot(2,3,2)
plot(ts,[c_ss,c(1:Tplot)],'Color',colors_list{2},'LineWidth',2)
xlabel('Year','Interpreter','Latex')
ylabel('Consumption','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')

subplot(2,3,5)
plot(ts,(1+[spread,s(1:Tplot)])-1,'Color',colors_list{2},'LineWidth',2)
xlabel('Year','Interpreter','Latex')
ylabel('Spread','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')

f.Position = [200 200 700 400]
saveas(gcf, "../OutputForPaper/model_sigma75.png")