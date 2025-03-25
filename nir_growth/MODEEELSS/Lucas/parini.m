% initial Parameter value
A=1;
beta=0.3;
gamma=0.3;
delta=0.1;
rho=0.05;
sigma=1.5;

% steady state growth rates
gr1=(1-beta)*(rho-delta)/(gamma-sigma*(1-beta+gamma));
gr2=gr1*((1-beta+gamma)/(1-beta));
