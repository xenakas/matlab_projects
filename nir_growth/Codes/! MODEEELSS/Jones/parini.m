% initial Parameter values
sigA=0.6;
sigL=0.6;
sigK=0.4;
delta=0.05;     
nPop=0.015;
etaA=0.6;      
etaL=0.5;      
etaLp=0.6;
rho=0.04;      
gammap=1;        
alphaF=1;
alphaJ=1;

% Calculation of steady state growth rates (they are not necessary in other
% models)
betaK=((1-etaA)*sigL+etaL*sigA)/((1-etaA)*(1-sigK));
betaA=etaL/(1-etaA);
