%relaxsetting.m

% number of initial boundary conditions ( = number of state variables)
n1=2;

%number of differential equations (Note: number of final boundary
%conditions equals n-n1)
n=4;

% number of static equations to be solved simultanously
n3=0;

%number of mesh points
M=40;

%number of variables altogether
N=n+n3;         %(Do not change)

%Normalization of variables
normal=[];

% calculation of steady state values
kss=1;
uss=1-gr1/delta;
css=(delta*(beta-gamma)/beta*uss+delta*(1-beta+gamma)/beta)*kss;
hss=(A*beta*kss^(beta-1)*uss^(1-beta)/(rho+gr2*sigma))^(1/(-1+beta-gamma));
 
%guess of final steady state values
y=ones(N,1);
y(1)=kss;    %k
y(2)=hss;    %h
y(3)=uss;    %u
y(4)=css;    %c

%Introduction of a shock (here, a deviation from the steady state
%relation in (k,h) induces transitional dynamics)
statev=ones(N,1);
statev(1)=1;
statev(2)=1.5;      

%--------------------------------------------------------------------------
%Specification, which differential equations are used for constructing the
%final boundary conditions. 
Endcond=[1 2];

tol=10^-9;       %tolerance for the Newton procedure
maxit=100;       %Maximum number of iterations
nu=0.2;          %Parameter for time transformation
damp=1;          %Dampening factor of the Newton procedure. The dampening factor will be 
dampfac=2;       %multiplied by the factor dampfac in every iteration until it equals 1
