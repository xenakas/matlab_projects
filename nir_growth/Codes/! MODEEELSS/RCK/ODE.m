%input of system of ordinary differential equations (as a vector)
 
funcODE=[ (k^alpha-c-(nPop+delta+gx)*k); ...         %ODE: k
          (c*(alpha*k^(alpha-1)-(delta+rho+gx*theta))/theta);  %ODE: c
             ];