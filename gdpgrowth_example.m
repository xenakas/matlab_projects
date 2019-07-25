%%%%% GDP growth as a Markov chain, example

P = [0.5,0.5,0.0;
     0.03,0.9,0.07;
     0.0,0.2,0.8];          %% n-by-n transition matrix
 
pi0 = [0,1,0];              %% initial probability distribution

mu  = 0.02;                 
sig = 0.04;

x   = [mu-sig;mu;mu+sig];   %% n-by-1 state

T   = 100;                  %% length of simulation

%%%%% Simulation

[chain,states] = simulate_markov(x,P,pi0,T);

logy = cumsum(chain);
ts   = cumsum(ones(T,1));


figure(1)
subplot(1,2,1)
plot(ts,chain,'bo-',ts,(mu-sig)*ones(T,1),'k--',ts,mu*ones(T,1),'k--',ts,(mu+sig)*ones(T,1),'k--')
title('GDP growth, a Markov chain')
subplot(1,2,2)
plot(logy,'b-')
title('log GDP')

