function y = errss_alt(par,P_ss,c_ss,M_ss,B_ss,r_ss,rb_ss)

global nu psi sigma 

alpha = par(1);
eta = par(2);

spread = (r_ss-rb_ss)/(1+r_ss);

iota = (alpha * (r_ss/(1+r_ss))^(1-nu)+...
    (1-alpha) * spread^(1-nu))^(1/(1-nu));

H = (iota*c_ss^(-sigma)/eta)^(-1/psi)*P_ss;

y(1) = M_ss - H*alpha *(r_ss/(1+r_ss)/iota)^(-nu);

y(2) = B_ss - H*(1-alpha)*(spread/iota)^(-nu);



% y = z(1)^2 + z(2)^2;