%relaxsetting.m

% number of initial boundary conditions ( = number of state variables)
n1=1;

%number of differential equations (Note: number of final boundary
%conditions equals n-n1)
n=2;

% number of static equations to be solved simultanously
n3=0;

%number of mesh points
M=60;

%number of variables altogether
N=n+n3;         %(Do not change)

%Normalization of variables
normal=[];  %no normalisation

% Steady state values
kss=((delta+rho+gx*theta)/alpha)^(1/(alpha-1));
css=kss^alpha-(nPop+delta+gx)*kss;

%guess of final steady state values
y=ones(N,1);
y(1)=kss;    %k
y(2)=css;    %c

%In case the shock consists of a reduction of state variables, enter this
%here
statev=ones(N,1);
statev(1)=0.1;      % The initial value of capital is at 10% of its steady state value

%--------------------------------------------------------------------------
%Specification, which differential equations are used for constructing the
%final boundary conditions. 
Endcond=0;

tol=10^-9;      %tolerance for the Newton procedure
maxit=50;       %Maximum number of iterations
nu=0.05;        %Parameter for time transformation
damp=1;         %Dampening factor of the Newton procedure. The dampening factor will be 
dampfac=2;      %multiplied by the factor dampfac in every iteration until it equals 1

