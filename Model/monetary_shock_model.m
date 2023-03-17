clear
close all
clc

global sigma nu psi

colors_list = {'#e41a1c','#377eb8', '#4daf4a','#984ea3','#ff7f00'};

%% Parameters

sigma = 5;
nu = 1;
psi = 1;                        % nu and psi are parameters for a general CES function for utility 
                                % from money and bonds. Set them to 1 for log-separable
delta_annual = 1/1.1;           % Discount factor
mubar_annual = 0.02;            % Steady state inflation
spread_annual = 0.08;           % Steady state spread
freq = 1;                       % Frequency of the model (1 for annual, 4 for quarterly, etc.)
delta = delta_annual^(1/freq);
mubar = (1+mubar_annual)^(1/freq)-1;
spread = (1+spread_annual)^(1/freq)-1;

%% Precision parameters
dist = 10;          
rL = 0.0001/freq;        % Low guess for initial interest rate
rH = 0.2/freq;           % High guess for initial interest rate
explode = rH;            % Highest value allowed during the path of interes rates
e = 1e-10;               % Convergence criterion
T = 150;                 % Horizon of simulation

%% Targets for steady state
B_ss = .5/.68;           % Target: public debt to consumption ratio
M_ss = .16/.68;          % Target: M1 to consumption ratio
c_ss_annual = 1;         % Normalization
P_ss = 1;                % Normalization
c_ss = c_ss_annual/freq;
r_ss = 1/delta*(1+mubar)-1;
rb_ss = r_ss - spread*(1+r_ss);

%% Solve for eta and alpha that match steady state targets

guesspar = [.25;.1];            

[par, e_par] = fsolve(@(e) err(e,P_ss,c_ss,M_ss,B_ss,r_ss,rb_ss),guesspar);
                               
alpha = par(1);
eta = par(2);       % alpga and eta are parameters of the CES utility of money and bonds


%% Path of shocks

M(1) = 1;           % On impact, M does not change
M(2) = 0.9836;      % M falls one period later
B(1) = 1.00914;     % Bonds rise on impact

M(3:T) = M(2);      % M stays constant after period 2

P(1) = P_ss;
B(1) = B(1)*B_ss;   % Rescale by steady state bond supply


M = M*M_ss.*(1+mubar).^(0:T-1);     % M scales up by steady state inflation

% Compute initial total liquidity H
H(1) = (alpha^(1/nu) * M(1)^((nu-1)/nu) + (1-alpha)^(1/nu) * B(1)^((nu-1)/nu))^(nu/(nu-1));


%% Solve forward. Guess a value for the zero-beta rate in period 1, solve forward and update

while dist > e
    % Update guess of initial r0
    r(1) = (rL +rH)/2;      

    % Compute cost of liquidity using quantities M, H and r0
    iota(1) = r(1)/(1+r(1))*alpha^(-1/nu)*(M(1)/H(1))^(1/nu);

    % Compute the spread using bond demand
    s(1) = iota(1)*(B(1)/H(1))^(-1/nu)*(1-alpha)^(1/nu);

    % Compute safe rate from r0 and the spread
    rb(1) = r(1) - s(1)*(1+r(1));
    
    % Compute consumption using liquidity demand
    c(1) = eta^(-1/sigma)*(H(1)/P(1))^(psi/sigma)*(iota(1))^(1/sigma);
    
    for t = 2:T
        % Bond supply to restore M/B to steady state
        B(t) = M(t)*B_ss/M_ss;

        % Price level fixed (follows steady state inflation)
        P(t) = P_ss*(1+mubar)^(t-1);

        % Compute total liquidity supply
        H(t) = (alpha^(1/nu) * M(t)^((nu-1)/nu) + (1-alpha)^(1/nu) * B(t)^((nu-1)/nu))^(nu/(nu-1));

        % Compute r0_{t+1} from difference equation that results from
        % replacing money demand into Euler equation
        r(t) = delta*r(t-1)/...
            ( (H(t)/H(t-1))^(psi-1/nu)*(M(t)/M(t-1))^(1/nu)*(P(t)/P(t-1))^(1-psi) - delta*r(t-1));	

        % Compute cost of liquidity using quantities M, H and r0
        iota(t) = r(t)/(1+r(t))*alpha^(-1/nu)*(M(t)/H(t))^(1/nu);

        % Compute the spread using bond demand
        s(t) = iota(t)*(B(t)/H(t))^(-1/nu)*(1-alpha)^(1/nu);

        % Compute safe rate from r0 and the spread
        rb(t) = r(t) - s(t)*(1+r(t));

        % Compute consumption using liquidity demand
        c(t) = eta^(-1/sigma)*(H(t)/P(t))^(psi/sigma)*(iota(t))^(1/sigma);
    end
    
    if real(r(T))>r_ss
        rH = r(1);              % If r0(T)>r_ss, initial guess was too high
    elseif max(real(r)) > explode
        rH = r(1);              % If r0(t) becomes too high as some point, initial guess was too high
    else
        rL = r(1);              % If r0(T)<r_ss, initial guess was too low
    end
       
    % Compute distance to convergence
    if max(real(r)) > explode
        dist = 10;              
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
dr0 = ((1+r(1))^freq - (1+r_ss)^freq)*100

% Safe rate
drb = ((1+rb(1))^freq - (1+rb_ss)^freq)*100

% Consumption level
dlogc = (log(c(1))-log(c_ss))*100 

% Consumption growth
dgc = (c(2)/c(1)-1)*100 


%% 3-period plots

Tplot = 3;

ts = linspace(0,Tplot,Tplot+1);

f = figure(1)

subplot(2,3,4)
plot(ts,[M_ss/(1+mubar),M(1:Tplot)]./(1+mubar).^(ts-1)/M_ss,'Color',colors_list{2},'LineWidth',2)
xlabel('Year','Interpreter','Latex')
ylabel('Money Supply','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')


subplot(2,3,1)
plot(ts,[B_ss/(1+mubar),B(1:Tplot)]./(1+mubar).^(ts-1)/B_ss,'Color',colors_list{2},'LineWidth',2)
xlabel('Year','Interpreter','Latex')
ylabel('Bond Supply','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')


subplot(2,3,3)
plot(ts,(1+[(1+r_ss)/(1+mubar)-1,r_real(1:Tplot)]).^freq-1,'Color',colors_list{2},'LineWidth',2)
xlabel('Year','Interpreter','Latex')
ylabel('Real Zero-Beta Rate','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')


subplot(2,3,6)
plot(ts,(1+[(1+rb_ss)/(1+mubar)-1,rb_real(1:Tplot)]).^freq-1,'Color',colors_list{2},'LineWidth',2)
xlabel('Year','Interpreter','Latex')
ylabel('Real Safe Rate','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')


subplot(2,3,2)
plot(ts,[c_ss,c(1:Tplot)],'Color',colors_list{2},'LineWidth',2)
xlabel('Year','Interpreter','Latex')
ylabel('Consumption','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')

subplot(2,3,5)
plot(ts,(1+[spread,s(1:Tplot)]).^freq-1,'Color',colors_list{2},'LineWidth',2)
xlabel('Year','Interpreter','Latex')
ylabel('Spread','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')

f.Position = [200 200 700 400]
saveas(gcf, "../Output/model.png")
