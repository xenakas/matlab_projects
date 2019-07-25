%%%%% Markov chain example

P = [0.7,0.2,0.1;
     0.0,0.5,0.5;
     0.0,0.9,0.1]; %% n-by-n transition matrix
 
pi0 = [1,0,0];     %% initial probability distribution

x   = [1;2;3];     %% n-by-1 state

T   = 20;          %% length of simulation

%%%%% Stationary distributions

%% By brute force....

pibars      = P^1000

pibar_brute = pibars(1,:)

%% By eigenvector decomposition

[V,D] = eig(P');
d     = diag(D);
[smallest,index] = min(abs(d-1)); %% find a unit eigenvalue
v     = V(:,index)

pibar_nice = (v/sum(v))'


%%%%% Simulation

[chain,states] = simulate_markov(x,P,pi0,T);

figure(1)
plot(chain,'bo-')
title('Markov chain')

