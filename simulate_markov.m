function [chain,state] = simulate_markov(x,P,pi0,T);

%%  x     = the quantity corresponding to each state, typical element x(i)
%%  P     = Markov transition matrix, typical element p(i,j) i,j=1,...n
%%  pi0   = probability distribution over initial state
%%  T     = number of periods to simulate
%%  
%%  chain = sequence of realizations from the simulation
%%  Modification of progam by L&S.


n = length(x); %% what is the size of the state vector?
E = rand(T,1); %% T-vector of draws from independent uniform [0,1]  

cumsumP = P*triu(ones(size(P)));
               %% creates a matrix whose rows are the cumulative sums of
               %% the rows of P               
               
%%%%% SET INITIAL STATE USING pi0               

E0   = rand(1,1);
ppi0 = [0,cumsum(pi0)];
s0   = ((E0<=ppi0(2:n+1)).*(E0>ppi0(1:n)))';
s    = s0; 

%%%%% ITERATE ON THE CHAIN

for t=1:T,
    state(:,t) = s;
    ppi        = [0,s'*cumsumP];
    s          = ((E(t)<=ppi(2:n+1)).*(E(t)>ppi(1:n)))';
end

chain = x'*state;    


