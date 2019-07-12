% MATLAB program solving a neoclassical stochastic 
% growth model by using the PEA with the moving bounds 
% as described in the article "Parameterized Expectations 
% Algorithm and the Moving Bounds" by Lilia Maliar and 
% Serguei Maliar, JBES, Volume 21, Nº 1, January 2003. 
% 
% The program includes two files: 
%
% 1. PEAmbound.m 
% 2. objective.m 

%----------------------------------------------------

clear all;

% 1. Initialize the model parameters
% ----------------------------------
alpha   = 0.33;         % Capital share 
delta   = .95;          % Discount factor
gam     = 3.0;          % Risk aversion parameter
d       = .02;          % Depreciation rate
sigma   = .01;          % Standard deviation for log noise
rho     = 0.95;         % Persistence of log technology shock
T       = 1000;         % Length of simulation

% 2. Allocate memory for the simulates series
% -------------------------------------------
tet  = zeros(T,1);      % Technology shocks
k    = zeros(T+1,1);    % Capital
c    = zeros(T,1);      % Consumption
e    = zeros(T-1,1);    % Conditional expectation

% 3. The steady state of k, c and the expectation
% -----------------------------------------------
ks = ( (1-delta+delta*d) / (alpha*delta) )^(1/(alpha-1) );
cs = ks^alpha - d*ks;												  
es = cs^(-gam)*( 1-d+alpha*ks^(alpha-1) );

% 4. Initial values of capital and technology shock
% -------------------------------------------------
k(1)   = ks;            % Initial value of capital 
tet(1) = 1;             % Initial value of shock

% 5. Draw the Technology Shocks
% -----------------------------
epsi = randn(T,1)*sigma;
for t = 2:T; tet(t) = tet(t-1)^rho*exp(epsi(t)); end;

% 6. Initialize the Algorithm Parameters
% --------------------------------------
beta    	= [log( es ); 0.00001; 0.00001];  	 % Initial polinomial coefficients  
crate    = 0.007;                             % Speed of moving the bound
criter  	= 1e-5;            	            	 % Convergence criterion
update  	= 1.0;            	                % Updating parameter for homotopy 

% 7. The Main Loop 
% -----------------             
iteration  = 0;                               % Initially, iteration is 0
dif	     = 2e-5;					             % Initally, criterion is not satified

while (dif > criter)|(hit==1)
up_bound  = ks*(2-exp(-crate*iteration));     % Upper bound
low_bound = ks*exp(-crate*iteration);         % Lower bound
hit       = 0;                                % Indicator, 1 (0) = bound is (not) hit ;  

% 7.1 Given 'beta', compute the time series
% -----------------------------------------

for t = 1:T
      uprime = exp( beta(1) + beta(2)*log(k(t)) + beta(3)*log(tet(t)));
      c(t) = ( delta*uprime )^(-1/gam);
      k(t+1) = k(t)^alpha*tet(t) - c(t) + (1-d)*k(t);
      if k(t+1) > up_bound  
          k(t+1) = up_bound; hit=1; 
      elseif k(t+1) < low_bound
          k(t+1) = low_bound; hit=1;
      end;
      c(t) = k(t)^alpha*tet(t)  + (1-d)*k(t)-k(t+1);
   end;
   
% 7.2 Given simulated time series, compute the expectation part
% -------------------------------------------------------------
   for t = 1:T-1
       e(t) = c(t+1)^(-gam)*(1-d+k(t+1)^(alpha-1)*alpha*tet(t+1));
   end;
    
% 7.3 Recompute 'beta' by using NLLS regression
% ---------------------------------------------
   x = [ones(T-1,1) log( k(1:T-1) ) log( tet(1:T-1) )];  % Regressors 
   ksi = nlinfit(x,e,'objective',beta);                  % NLLS regression
   iteration                                             % Display iteration
   dif = norm(beta-ksi)                                  % Display difference between 
   beta = update*ksi + (1-update)*beta;                  % Update the coefficients (homotopy)
   iteration = iteration+1;			                     % Next iteration
      
end;


% 8. Plot the time series solution y, c and k 
% -------------------------------------------
time=(1:1:T);                         
subplot(3,1,1);
plot (time,k(1:T,1)), xlabel('t'), ylabel('Capital')
title('Time series solution');
subplot(3,1,2);
plot (time,c), xlabel('t'), ylabel('Consumption')
subplot(3,1,3);
plot (time,tet.*k(1:T,1).^alpha), xlabel('t'), ylabel('Output')