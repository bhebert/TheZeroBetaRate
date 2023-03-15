function y = err(par,P_ss,c_ss,M_ss,B_ss,r_ss,rb_ss)

% This function gives the error in the money and bond demand equations
% given parameters, money and bond supply zero-beta rate and safe rate

global nu psi sigma 

alpha = par(1);
eta = par(2);

% Compute spread
spread = (r_ss-rb_ss)/(1+r_ss);

% Compute total cost of liquidity
iota = (alpha * (r_ss/(1+r_ss))^(1-nu)+...
    (1-alpha) * spread^(1-nu))^(1/(1-nu));

% Compute composite money supply
H = (iota*c_ss^(-sigma)/eta)^(-1/psi)*P_ss;

% Error in money demand equation
y(1) = M_ss - H*alpha *(r_ss/(1+r_ss)/iota)^(-nu);

% Error in bond demand equation
y(2) = B_ss - H*(1-alpha)*(spread/iota)^(-nu);
