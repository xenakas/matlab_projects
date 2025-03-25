function funcODE = funcODE(t,x)
% System file of the relaxation algorithm by Trimborn, Koch, and Steger.
% Copyright by Trimborn, Koch, Steger, 2008

global alpha delta rho theta    lambda  determ   gamma2  c_tilde_spline  direction  c_tilde_spline 


%list of variables to extract
k=x(1,:);
c=x(2,:);

% solves either upper arm (direction=1) or lower arm (direction=0) of the policy function
if direction==1
    k=1./k;
    c=1./c;    
end

% factor rewards
r=alpha*k.^(alpha-1);
w=(1-alpha)*k.^alpha;

% specifies the deterministic system
if determ==1
    funcODE=[ (k^alpha-delta*k-c); ...          %ODE: k
              (c.*(r-delta-rho)/theta);               %ODE: c
             ];
% specifies the conditional deterministic system (conditioned on no stochastics)
elseif determ==2  
    funcODE=[ (k^alpha-delta*k-c); ...                                             %ODE: k
              (c.*(r-delta-rho-lambda+lambda*(1-gamma2)*(ppval(c_tilde_spline,k))^(-theta))/theta);  %ODE: c
             ];
end

% transforms ODEs for computation of upper arm (transformation) or lower arm (no transformation) of the policy function
if direction==1
    funcODE(1)=-1/k^2*funcODE(1);
    funcODE(2)=-1/c^2*funcODE(2);
end

% only applicable for backward integration
funcODE=funcODE*-1;
