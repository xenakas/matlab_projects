%relaxsetting.m

% number of initial boundary conditions ( = number of state variables)
n1=1;

%number of differential equations (Note: number of final boundary
%conditions equals n-n1)
n=1;

% number of static equations to be solved simultanously
n3=2;       % To add static equations is optional

%number of mesh points
M=50;

%number of variables altogether
N=n+n3;         %(Do not change)

%Normalization of variables
normal=[];

% calculation of steady state values
kss=(s/(gx+nPop+delta))^(1/(1-a));
css=(1-s)*kss^a;
yss=kss^a;

%guess of final steady state values
y=ones(N,1);
y(1)=kss;   %k
y(2)=css;   %c
y(3)=yss;   %y

%In case the shock consists of a reduction of state variables, enter this
%here
statev=0;       %Use initbound.m to determine initial boundary conditions

%--------------------------------------------------------------------------
%Specification, which differential equations are used for constructing the
%final boundary conditions. 
Endcond=[0];

tol=10^-9;      %tolerance for the Newton procedure
maxit=50;       %Maximum number of iterations
nu=0.04;        %Parameter for time transformation
damp=0.3;       %Dampening factor of the Newton procedure. The dampening factor will be 
dampfac=2;      %multiplied by the factor dampfac in every iteration until it equals 1
