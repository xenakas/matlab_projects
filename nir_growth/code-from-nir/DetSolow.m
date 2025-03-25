%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	DETERMINISTIC GROWTH MODEL SOLVED WITH DYNAMIC PROGRAMMING 
%
%
%   max \sum_t beta^t ln c_t
%
%   subject to c_t + k_{t+1} = exp(sigma) k_t^{alpha} 
%
%
%  George Hall
%  Brandeis University
%  Econ 303: Advanced Macroeconomics I
%  Fall 2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
%! rm det_dp.out
diary det_dp.out; 
disp('DETERMINISTIC GROWTH MODEL SOLVED WITH DP');
disp('');
%
%  set parameter values
%
alpha  = 0.30;            % capital's share of income
beta   = 1.03^(-0.25);    % subjective discount factor 
sigma  = 0.22;            % TFP = exp(sigma)              
%
%   form capital grid
%   
maxkap = 1.00;                  % maximum value of capital grid   
inckap = 0.001;               % size of capital grid incremenats
minkap = 0;                  % minimum value of capital grid   
nkap   = round((maxkap-minkap)/inckap+1); % number of grid points
kgrid = [ minkap:inckap:maxkap ]';
%
%  compute the true value function from the analytical solution
%
true_val_func = 1/(1-beta)*(log(exp(sigma)*(1-beta*alpha)) + beta*alpha/(1-beta*alpha)*log(exp(sigma)*beta*alpha)) + alpha/(1-alpha*beta)*log(kgrid); 
%
%  compute true decision rule from analytical solution
%
true_dec = alpha*beta*exp(sigma)*kgrid.^(alpha);
%
%  compute the steady state of the capital stock from analytical solution
%
ssk = (alpha*beta*exp(sigma))^(inv(1-alpha));
%
%
%  Now compute the solution numerically ...
% 
%  tabulate the utility function such that for zero or negative
%  consumption utility remains a large negative number so that
%  such values will never be chosen as utility maximizing      
%

kapp = repmat(kgrid,1,nkap);
kap = repmat(kgrid,1,nkap)';
%
%  impose budget constraint
%
cons = exp(sigma)*kap.^alpha - kapp;

clear kapp kap
%
%  rule out negative consumption
%
cons(find(cons<=0)) = NaN;
%
%  Assume log utility 
%
util =  log(cons);
%
%  if consumption is negative, set utility to negative infinity
%
util(find(isnan(util))) = -inf;

clear cons

%
%  initialize some variables
%
iter    = 0;
v       = -76*ones(nkap,1);
decis   = zeros(nkap,1);
test    = 10;
[rs,cs] = size(util);
%
%
%  iterate on Bellman's equation and get the decision 
%  rules and the value function at the optimum         
%
format short g
while test > 1e-7;

    [tv,tdecis]=max(util + beta*repmat(v,1,nkap));

    tv=tv';

    test=max(abs((tv-v)./v));
    iter = iter + 1;
    %
    %  plot the value function as it is converging
    %
    if ( (iter < 10) | (iter/10) == floor(iter/10)); 
        disp([ iter test ]); 
        figure(1)
        plot(kgrid,v,kgrid,true_val_func)
        ititle = ['value function at iteration ' num2str(iter) ];
        title(ititle);
        axis([ 0 1 -79 -75 ]);
        xlabel('capital stock');
        ylabel('value function');

        if(iter < 100); 
            pause(1)
        else
            pause(.05)
        end
    end;

    v=tv;
    decis=tdecis';
end;
%
%  decis is decision rule computed numerically by the computer
%
decis=(decis-1)*inckap + minkap;
%
%   plot the decision rule and value function
%

figure(2)
plot(kgrid,decis,'b-',kgrid,true_dec,'r:',kgrid,kgrid,'g',ssk,ssk,'*')
title('decision rules - deterministic dynamic programming');
xlabel('capital this period');
ylabel('capital next period');
print kdecis.ps

figure(3)
plot(kgrid,v)
title('converged value function - deterministic dynamic programming');
xlabel('capital stock');
ylabel('value function');
print value.ps


diary off;