function [t, x] = relax(funcODE, funcSTAT, funcINI, funcfinal, n, n1, n3, nu, y, M, start, Endcond, maxit, tol, damp, dampfac)                           
% Version 3.1
% RELAXATION algorithm to solve infinite-horizon continuous time models.
%
% Description of procedure:
% Trimborn, Koch, Steger (2007) Multidimensional Transitional Dynamics: A
% Simple Numerical Procedure, forthcoming in Macroeconomic Dynamics
% 
% Copyright by Trimborn, Koch, Steger, 2008
%
% For further information contact Timo Trimborn, University of Hamburg
% or visit http://www.relaxation.uni.siegen.de
%
% COMMAND:
% [t, x] = relax(@funcODE, @funcSTAT, @funcINI, @funcfinal, n, n1, n3, nu, y, M, start, Endcond, maxit, tol, damp, dampfac)
% 
% INPUT ARGUMENTS
% funcODE: function which returns the right hand side of the 
%     differential equations whereas time and vector of variales are inputs
% funcSTAT: function which returns the residual of the static equations
%     whereas time and vector of variables are inputs
% funcINI: function which returns the residual of the initial boundary conditions
%     (In case you want to determine them otherwise, nevertheless input @funcINI)
% funcfinal: function which returns the residual of the final boundary conditions
%     (In case you want to determine them otherwise, nevertheless input @funcfinal)
% n: total number of differential equations
% n1: number of initial boundary conditions
% n3: (optional, default value is 0) number of static equations
% nu: (optional, default value is 0.1) parameter for time transformation
% y: (optional, default value is vector 1) vector of dimension M*Nx1, guess of the 
%     time path of the variables. the first N variables are the guess at time 0
%     and so on whereas the time vector is determined through M and nu     
% M: (optional, default value is 50) number of mesh points
% start: (optional, default is start=0) vector with the values of the state variables 
%     at time 0. If input is start=0, external file funcINI will be used to determine 
%     the initial boundary conditions
% Endcond: (optional, default is Endcond=0) vector with determines, which differential
%     equations will be used for constructing the final boundary conditions. If input
%     is Endcond=0, external file funcfinal will be used to determine the final 
%     boundary conditions
% maxit: (optional, default is maxit=50) maximum number of iterations
% tol: (optional, default is tol=10^-9) tolerance for Newton-procedure
% damp: (optional, default is damp=1) Dampening factor of the Newton procedure. The
%     dampening factor will be multiplied by the factor dampfac in every iteration
%     until it equals 1. Therefore damp=1 means no dampening.
% dampfac: (optional, default is dampfac=2)
%
% OUTPUT ARGUMENTS:
% t: time vector
% x: solution vector of variables through time: column i are the values at time t(i)    

%**************************************************************************
%--------------------------------------------------------------------------
N=n+n3;                    %calculation of total dimension
if nargin < 16, dampfac=2; end                        %default values
if nargin < 15, damp=1; end
if nargin < 14, tol=10^-9; end
if nargin < 13, maxit=50; end
if nargin < 12, Endcond=0; end
if nargin < 11, start=0; end
if nargin < 10, M=50; end
if nargin < 9, y=ones(N*M,1); end
if nargin < 8, nu=0.1; end
if nargin < 7, n3=0; end
crit=(y(end-N+1:end).*y(end-N+1:end)<eps*10)*10^-5;   %The convergence criteria will be (deltay'.*deltay) normalized by the final steady state vector,
crit=y(end-N+1:end)+crit;                             %therefore for components which are zero an increment has to be added to the final vector 
it=0;                       %reset number of iterations
E=zeros(N*M,1);             %This vector contains the residuals from the nonlinear equations
fac=[];fac2=[];fac3=[];     %This is a storing vector for the numerical calculation of the jacobian
Jac=zeros(N*(M-1),N);       %Later on, in this vector the sections of the jacobian matrix will be stored
KonvKrit=tol+1;
%start of iteration+++++++++++++++++++++++++++++++++++++++++++++++
while KonvKrit>tol & it<maxit    
    if it==0; disp(['Start of main loop:']); end; 
    it=it+1;    
    
    %Input of initial conditions and first static equation in E
    if start==0
        E(1:n1)=feval(funcINI,0,y(1:N));
    else
        E(1:n1)=y(1:n1)-start(1:n1);
    end
    if n3>0
        E(n1+1:n1+n3)=feval(funcSTAT,0,y(1:N));    
        %creation of the first part of the jacobian due to the initial
        %conditions and the first static equation
        [Jacstat,fac2]=numjac(funcSTAT,0,y(1:N),feval(funcSTAT,0,y(1:N)),sqrt(eps)*ones(N,1),fac2,0);
    else
        Jacstat=[];
    end
    if start==0
        [Jacini,fac3]=numjac(funcINI,0,y(1:N),feval(funcINI,0,y(1:N)),sqrt(eps)*ones(N,1),fac3,0);
        Jacini=Jacini(1:n1,:);
    else
        Jacini=[eye(n1) zeros(n1,N-n1)];
    end
    Jacsml=[Jacini ; Jacstat];
    for i=1:M-1
        %numerical evaluation of H, SL and SR
        E(n1+n3+(i-1)*N+1:n1+n3+i*N)=feval(@funcH,i/(M-1),y(i*N+1:i*N+N),y((i-1)*N+1:(i-1)*N+N),(i-1)/(M-1),n,N,nu,funcODE,funcSTAT);    %function H is evaluated
        [SRnum,fac]=numjac(@funcH,i/(M-1),y(i*N+1:i*N+N),feval(@funcH,i/(M-1),y(i*N+1:i*N+N),y((i-1)*N+1:(i-1)*N+N),(i-1)/(M-1),n,N,nu,funcODE,funcSTAT),...
            sqrt(eps)*ones(N,1),fac,0,[],[],y((i-1)*N+1:(i-1)*N+N),(i-1)/(M-1),n,N,nu,funcODE,funcSTAT);    %Jacobian block SR is calculated
        SLnum=SRnum-2*eye(size(SRnum));                     %Jacobian block SL can be computed out of SR
        SLnum(N-n3+1:N,:)=0;
        Jactmp=[Jacsml zeros(size(Jacsml)); SLnum SRnum];   %The block which has to be modified is composed 
        [Jactmp E((i-1)*N+1:i*N+n1+n3)]=feval(@rrefmod,Jactmp , E((i-1)*N+1:i*N+n1+n3));   %Transformation of actual part of the jacobian
        Jac((i-1)*N+1:i*N,:)=Jactmp(1:N , N+1:2*N);         %relevant areas are stored
        Jacsml=Jactmp(N+1:N+n1+n3,N+1:2*N);                 %One part is stored for the next passage of the loop
    end
    
    %input of final boundary conditions in E
    End=feval(funcODE,inf,y(N*M-N+1:N*M));
    if n>n1 & Endcond~=0
        E(N*M-n+n1+1:N*M)=End(Endcond);  
    elseif n>n1
        E(N*M-n+n1+1:N*M)=feval(funcfinal,inf,y(N*M-N+1:N*M));
    end
              
    %Calculation of part of the jacobian due to final boundary conditions
    if Endcond~=0 | n==n1
        [Jacend,fac3]=numjac(funcODE,inf,y(N*(M-1)+1:N*M),feval(funcODE,inf,y(N*(M-1)+1:N*M)),sqrt(eps)*ones(N,1),fac3,0);
    else
        [Jacend,fac3]=numjac(funcfinal,inf,y(N*(M-1)+1:N*M),feval(funcfinal,inf,y(N*(M-1)+1:N*M)),sqrt(eps)*ones(N,1),fac3,0);
    end
    
    %Transforming final part of the jacobian
    if n>n1 & Endcond~=0
        Jactmp=[Jacsml ; Jacend(Endcond,:)];
    elseif n>n1
        Jactmp=[Jacsml ; Jacend];
    else
        Jactmp=Jacsml;
    end
    [Jactmp E(N*(M-1)+1:N*M)]=feval(@rrefmod,Jactmp , E(N*(M-1)+1:N*M));

    %Calculation of correction vector and correction of last block
    deltay=E(N*M-N+1:N*M);
    y(N*M-N+1:N*M) = y(N*M-N+1:N*M) - damp*deltay;
    %Calculation of convergence criteria
    KonvKrit=sqrt((deltay./(crit))'*(deltay./(crit)));
    %Calculation of correction vector and correction of all other blocks
    for i=M-1:-1:1
        deltay=E((i-1)*N+1:i*N)-Jac((i-1)*N+1:i*N,:)*deltay;
        y((i-1)*N+1:i*N) = y((i-1)*N+1:i*N) - damp*deltay;
        KonvKrit=KonvKrit+sqrt((deltay./(crit))'*(deltay./(crit))); %convergence criteria
    end        
    %Calculation of convergence criteria
    KonvKrit=KonvKrit/(M*N);
    if damp<1
        disp(['Iteration number: ',num2str(it),'  convergence criteria: ',num2str(KonvKrit),'  dampening factor: ',num2str(damp)]);
    else
        disp(['Iteration number: ',num2str(it),'  convergence criteria: ',num2str(KonvKrit)]);
    end
    damp=min(dampfac*damp,1);
    y=real(y);
end;

if KonvKrit > tol | not(isfinite(KonvKrit))        %Test of convergence
    disp(['No convergence!']);
elseif KonvKrit~=0
    disp(['Convergence achieved.']);disp([' ']);
end
%End of iteration++++++++++++++++++++++++++++++++++++++++++++++++++++

%Storing of variables in x and creation of actual time t
x=zeros(N,M);
t=zeros(1,M);
tau=[0:1/(M-1):1];
for i=1:M
    x(:,i)=y((i-1)*N+1:i*N);
    if tau(i)~=1
        t(i)=tau(i)/nu/(1-tau(i));
    else
        t(i)=inf;
    end;
end;

%**************************************************************************
%--------------------------------------------------------------------------

function funcH = funcH(taua,ya,yb,taub,n,N,nu,funcODE,funcSTAT)
%ya: variables at knot k+1, yb: variables at knot k
%taua: time at k+1, taub: time at k
%The residuals of the differential equations and the static equations are composed to one
%vector. The differential equation is tranformed into a difference equation
%(midpoint rule)

time=((taua+taub)/2);
if taua<1
    funcH=[ya(1:n)-yb(1:n)-(taua-taub)*feval(funcODE,time/(nu*(1-time)),(ya+yb)/2)*1/nu/(1-time)^2; ...
       feval(funcSTAT,taua/(nu*(1-taua)),ya);];
else
    funcH=[ya(1:n)-yb(1:n)-(taua-taub)*feval(funcODE,time/(nu*(1-time)),(ya+yb)/2)*1/nu/(1-time)^2; ...
       feval(funcSTAT,inf,ya);];
end