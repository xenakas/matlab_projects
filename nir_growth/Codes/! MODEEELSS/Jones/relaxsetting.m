%relaxsetting.m

% number of initial boundary conditions ( = number of state variables)
n1=2;

%number of differential equations (Note: number of final boundary
%conditions equals n-n1)
n=4;

% number of static equations to be solved simultanously
n3=1;

%number of mesh points
M=100;

%number of variables altogether
N=n+n3;         %(Do not change)

%Normalization of variables
normal=[];

%guess of final steady state values
y=ones(N,1);
y(1)=600;    %k
y(2)=900;    %a
y(3)=600;    %c
y(4)=3;      %v
y(5)=0.9;    %theta

%In case the shock consists of a reduction of state variables, enter this
%here
statev=ones(N,1);       %Here a specific shock is simulated, specified in shock.m 

%--------------------------------------------------------------------------
%Specification, which differential equations are used for constructing the
%final boundary conditions. 
Endcond=[n1+1:n];

tol=10^-9;      %tolerance for the Newton procedure
maxit=50;       %Maximum number of iterations
nu=0.04;        %Parameter for time transformation
damp=1;         %Dampening factor of the Newton procedure. The dampening factor will be 
dampfac=2;      %multiplied by the factor dampfac in every iteration until it equals 1
