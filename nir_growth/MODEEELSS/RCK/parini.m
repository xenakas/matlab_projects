% initial Parameter values
alpha=0.3;
delta=0.05;
rho=0.02;
nPop=0.01;
gx=0.00;
theta =(delta+rho)/(alpha*(delta+nPop+gx)-gx); %Select theta so that the saving rate is constant
