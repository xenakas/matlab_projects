function findss=findss(k)

global   delta rho theta lambda  gamma2 alpha2 c_tilde_spline  c_tilde_spline 

r=alpha2*k^(alpha2-1)-delta;
findss=((r-rho-lambda+lambda*(1-gamma2)*(ppval(c_tilde_spline,k))^(-theta))/theta);

