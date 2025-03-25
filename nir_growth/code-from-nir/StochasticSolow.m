%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A Matlab program to solve a stochastic growth model via value function iteration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
format short g
! rm stoch_vfi.out
diary stoch_vfi.out; 
disp('A SIMPLE STOCHASTIC GROWTH MODEL');
disp('');
%
%  set parameter values
r_0 = 1.2; 

alpha  = 0.40;              % production parameter
beta   = 0.95;              % subjective discount factor 
prob   = [ .5 .5; .5 .5];   % prob(i,j) = probability (A(t+1)=Aj | A(t) = Ai)
delta  = .90;               % 1 - depreciation rate
A_high = 1.5;               % high value for technology
A_low  = 0.5;               % low value for technology
%
%   form capital grid
%   
mink =   0.01;                     % minimum value of the capital grid
maxk =  25.01;                     % maximum value of the capital grid   
ink  =  0.05;                     % size of capital grid increments
nk   = round((maxk-mink)/ink+1);   % number of grid points
kgrid = [ mink:ink:maxk ]';
% 
%  tabulate the utility function such that for zero or negative
%  consumption utility remains a large negative number so that
%  such values will never be chosen as utility maximizing      
%


kapp = repmat(kgrid,1,nk);
kap = repmat(kgrid,1,nk)';

cons1 = A_high*kap.^alpha  + delta*kap - kapp;
cons2 = A_low*kap.^alpha + delta*kap - kapp;

clear kapp kap

cons1(find(cons1<=0)) = NaN;
cons2(find(cons2<=0)) = NaN;

util1 =  log(cons1);
util2 =  log(cons2);

util1(find(isnan(util1))) = -inf;
util2(find(isnan(util2))) = -inf;

clear cons1 cons2
%
%  initialize some variables
%
v       = repmat(0,nk,2);
decis   = repmat(0,nk,2);
metric  = 10;
iter = 0;
tme = cputime;
[rs,cs] = size(util1);
%
%  iterate on Bellman's equation and get the decision 
%  rules and the value function at the optimum         
%
while metric > 1e-7;

  [tv1,tdecis1]=max(util1 + beta*repmat(v*prob(1,:)',1,nk));
  [tv2,tdecis2]=max(util2 + beta*repmat(v*prob(2,:)',1,nk));
  
  tdecis=[tdecis1' tdecis2'];
  tv=[tv1' tv2'];
  
  metric=max(max(abs(tv-v)./tv));
  v=tv;
  decis=tdecis;
iter = iter+1;
end;
disp('fixed point solved via value function iteration took');
disp([ iter ]);
disp('iterations and');
disp([ cputime-tme ]);
disp('seconds');

decis=(decis-1)*ink + mink;
%
%   form transition matrix
%   trans is the transition matrix from state at t (row)
%   to the state at t+1 (column) 
%   The eigenvector associated with the unit eigenvalue
%   of trans' is  the stationary distribution. 
% 
g2=sparse(cs,cs);
g1=sparse(cs,cs);
for i=1:cs
    g1(i,tdecis1(i))=1;
    g2(i,tdecis2(i))=1;
end
trans=[ prob(1,1)*g1 prob(1,2)*g1; prob(2,1)*g2 prob(2,2)*g2];
trans= trans';
probst = (1/(2*nk))*ones(2*nk,1);
test = 1;
while test > 10^(-8);
   probst1 = trans*probst;
   test=max(abs(probst1-probst));
   probst = probst1;
end;
%
%   vectorize the decision rule to be conformable with probst
%   calculate mean level of capital
%
kk=decis(:);
meanK=probst'*kk;
%
%  calculate measure over (k,A) pairs
%  lambda has same dimensions as decis
%
lambda=zeros(cs,2);
lambda(:)=probst;
%
%   calculate stationary distribution of capital 
%
probk=sum(lambda');     
probk=probk';
%
%   print out results
%
disp('PARAMETER VALUES');
disp('');
disp('    alpha      beta       '); 
disp([ alpha beta ]);
disp(''); 
disp('RESULTS ');
disp('');
disp('      mean of K ');
disp([ meanK ]);
%
%    simulate life histories of the agent
%
disp('SIMULATING LIFE HISTORY');
kgrid = [ (0:ink:maxk)' ];  % capital grid  
kmark = 10;
k = kgrid(kmark,1);        % initial level of assets
n = 100;                   % number of periods to simulate
s0 = 1;                    % initial state 
states   = zeros(n-1,2);
controls = zeros(n-1,2);
[chain,state] = markov(prob,n,s0);
for i = 1:n-1;
    if chain(i) == 1;
       kprime = decis(kmark,1);
       invest = kprime - delta*k; 
       cons   = A_high*k^(alpha) - invest; 
       kmark = tdecis(kmark,1);
    elseif chain(i) == 2;
       kprime = decis(kmark,2);
       invest = kprime - delta*k; 
       cons   = A_low*k^(alpha) - invest; 
       kmark = tdecis(kmark,2);
    else;
      disp('something is wrong with chain');                                                                                                                                                                             
    end;
    states(i,:) = [ k chain(i) ];
    controls(i,:) = [ cons kprime ];
    k = kprime;
end;


figure(1)
plot(kgrid',v(:,1),'-',kgrid',v(:,2),':');
title('STOCH GROWTH MODEL: VALUE FUNCTION');
print value.ps

figure(2)
plot(kgrid',decis(:,1),'.',kgrid',decis(:,2),':',kgrid',kgrid','-');
title('STOCH GROWTH MODEL: POLICY FUNCTION');
axis([ 0 maxk 0 maxk ]);
print policy.ps

figure(3)
plot((1:n-1)',controls(:,1));
title('STOCH GROWTH MODEL: CONSUMPTION');
print consum.ps

figure(4)
plot((1:n-1)',controls(:,2));
title('STOCH GROWTH MODEL: INVESTMENT');
print consum.ps

figure(5)
plot(kgrid,probk);
title('DISTRIBUTION OF CAPITAL');
xlabel('CAPITAL');
ylabel('FRACTION OF AGENTS');
print capdist.ps






kapp = repmat(kgrid,1,nk);
kap = repmat(kgrid,1,nk)';

cons1 = A_high*kap.^alpha  + delta*kap - kapp;
cons2 = A_low*kap.^alpha + delta*kap - kapp;

clear kapp kap

cons1(find(cons1<=0)) = NaN;
cons2(find(cons2<=0)) = NaN;

util1 =  log(cons1);
util2 =  log(cons2);

util1(find(isnan(util1))) = -inf;
util2(find(isnan(util2))) = -inf;

clear cons1 cons2
%
%  initialize some variables
%
v       = repmat(0,nk,2);
decis   = repmat(0,nk,2);
metric  = 10;
iter = 0;
tme = cputime;
[rs,cs] = size(util1);
%
%  iterate on Bellman's equation and get the decision 
%  rules and the value function at the optimum         
%
while metric > 1e-7;

  [tv1,tdecis1]=max(util1 + beta*repmat(v*prob(1,:)',1,nk));
  [tv2,tdecis2]=max(util2 + beta*repmat(v*prob(2,:)',1,nk));
  
  tdecis=[tdecis1' tdecis2'];
  tv=[tv1' tv2'];
  
  metric=max(max(abs(tv-v)./tv));
  v=tv;
  decis=tdecis;
iter = iter+1;
end;
disp('fixed point solved via value function iteration took');
disp([ iter ]);
disp('iterations and');
disp([ cputime-tme ]);
disp('seconds');

decis=(decis-1)*ink + mink;
%
%   form transition matrix
%   trans is the transition matrix from state at t (row)
%   to the state at t+1 (column) 
%   The eigenvector associated with the unit eigenvalue
%   of trans' is  the stationary distribution. 
% 
g2=sparse(cs,cs);
g1=sparse(cs,cs);
for i=1:cs
    g1(i,tdecis1(i))=1;
    g2(i,tdecis2(i))=1;
end
trans=[ prob(1,1)*g1 prob(1,2)*g1; prob(2,1)*g2 prob(2,2)*g2];
trans= trans';
probst = (1/(2*nk))*ones(2*nk,1);
test = 1;
while test > 10^(-8);
   probst1 = trans*probst;
   test=max(abs(probst1-probst));
   probst = probst1;
end;
%
%   vectorize the decision rule to be conformable with probst
%   calculate mean level of capital
%
kk=decis(:);
meanK=probst'*kk;
%
%  calculate measure over (k,A) pairs
%  lambda has same dimensions as decis
%
lambda=zeros(cs,2);
lambda(:)=probst;
%
%   calculate stationary distribution of capital 
%
probk=sum(lambda');     
probk=probk';
%
%   print out results
%
disp('PARAMETER VALUES');
disp('');
disp('    alpha      beta       '); 
disp([ alpha beta ]);
disp(''); 
disp('RESULTS ');
disp('');
disp('      mean of K ');
disp([ meanK ]);
%
%    simulate life histories of the agent
%
disp('SIMULATING LIFE HISTORY');
kgrid = [ (0:ink:maxk)' ];  % capital grid  
kmark = 10;
k = kgrid(kmark,1);        % initial level of assets
n = 100;                   % number of periods to simulate
s0 = 1;                    % initial state 
states   = zeros(n-1,2);
controls = zeros(n-1,2);
[chain,state] = markov(prob,n,s0);
for i = 1:n-1;
    if chain(i) == 1;
       kprime = decis(kmark,1);
       invest = kprime - delta*k; 
       cons   = A_high*k^(alpha) - invest; 
       kmark = tdecis(kmark,1);
    elseif chain(i) == 2;
       kprime = decis(kmark,2);
       invest = kprime - delta*k; 
       cons   = A_low*k^(alpha) - invest; 
       kmark = tdecis(kmark,2);
    else;
      disp('something is wrong with chain');                                                                                                                                                                             
    end;
    states(i,:) = [ k chain(i) ];
    controls(i,:) = [ cons kprime ];
    k = kprime;
end;


figure(1)
plot(kgrid',v(:,1),'-',kgrid',v(:,2),':');
title('STOCH GROWTH MODEL: VALUE FUNCTION');
print value.ps

figure(2)
plot(kgrid',decis(:,1),'.',kgrid',decis(:,2),':',kgrid',kgrid','-');
title('STOCH GROWTH MODEL: POLICY FUNCTION');
axis([ 0 maxk 0 maxk ]);
print policy.ps

figure(3)
plot((1:n-1)',controls(:,1));
title('STOCH GROWTH MODEL: CONSUMPTION');
print consum.ps

figure(4)
plot((1:n-1)',controls(:,2));
title('STOCH GROWTH MODEL: INVESTMENT');
print consum.ps

figure(5)
plot(kgrid,probk);
title('DISTRIBUTION OF CAPITAL');
xlabel('CAPITAL');
ylabel('FRACTION OF AGENTS');
print capdist.ps;


