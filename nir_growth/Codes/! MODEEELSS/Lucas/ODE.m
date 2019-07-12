%input of system of ordinary differential equation (as a vector)

funcODE=[      (A*k^(beta)*u^(1-beta)*h^(1-beta+gamma)-gr2*k-c); ...                     %ODE k;
               (h*(delta*(1-u)-gr1)); ...   %ODE: h            
               (delta*(beta-gamma)/beta*u^2+delta*(1-beta+gamma)/beta*u-c/k*u);                     %ODE:u  
               (c/sigma*(A*beta*k^(beta-1)*h^(1-beta+gamma)*u^(1-beta)-rho-gr2*sigma));  ...      %ODE: c
             ];
              
