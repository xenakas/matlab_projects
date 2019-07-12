%input of system of ordinary differential equations (as a vector)

funcODE=[ (alphaF*a^sigA*theta^sigL*k^sigK-c-delta*k-betaK*nPop*k); ...         %ODE: k
          (alphaJ*a^(etaA)*(1-theta)^etaL-betaA*nPop*a); ...                     %ODE a;
          (c/gammap*(sigK^2*alphaF*a^sigA*theta^sigL*k^sigK/k -delta-rho+gammap*nPop-nPop) - betaK*nPop*c); ...   %ODE: c
          (v*(sigK^2*alphaF*a^sigA*theta^sigL*k^sigK/k -delta-(betaK-betaA)*nPop) ...
          -(1-sigK)*sigK*alphaF*a^sigA*theta^sigL*k^sigK/a);              %ODE: v
             ];
             