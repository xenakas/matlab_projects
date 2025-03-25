function [guess, start, errorcode] = initrelax(funcODE, funcSTAT, n, n1, n3, nu, y, M, statev, t)
% Version 3.1
% Copyright by Trimborn, Koch, Steger, 2008
%
% COMMAND:
% [guess, start, errorcode] = initrelax(@funcODE, @funcSTAT, n, n1, n3, nu, y, M, statev, t)
% 
% INPUT ARGUMENTS:
% funcODE: function which returns the right hand side of the 
%     differential equations whereas time and vector of variales are inputs
% funcSTAT: function which returns the residual of the static equations
%     whereas time and vector of variables are inputs
% n: number of differential equations
% n1: number of initial boundary conditions
% n3: (optional, default value is 0) number of static equations
% nu: (optional, default value is 0.1) parameter for time transformation
% y: (optional, default value is vector 1) initial guess. If no time dependent guess is made, 
%     this should be a column vector, representing a guess for the steady state values of 
%     the variables (dimension n+n3). If a time dependent guess is made, this should be a
%     matrix with the same number of columns as t, such that each column i is a guess for the variables at time t(i)   
% M: (optional, default value is 50) number of mesh points
% statev: (optional, default value is 0) if statev=0, initial boundary conditions will be determined in the file initbound.m.
%     if statev=ones(N,1), the initial steady state will 
%     be determined via shock.m. If this is a vector of positiv numbers, the state 
%     variables will be multiplied by the correspondent entry
% t: (optional) if a time vector is supplied, y should match for a time dependent guess (see above)
%
% OUTPUT ARGUMENTS:
% guess: initial guess, suitable for relax.m
% start: vector which indicates treatment of initial boundary conditions
% errorcode: 0 no error 
%            1 not enough input arguments
%            2 funcODE or funcSTAT return non real values or vectors in wrong dimension
%            3 fsolve does not converge
%
% For further information contact Timo Trimborn, University of Hannover
% or visit http://www.relaxation.uni.siegen.de
    
N=n+n3;
errorcode=0;
if nargin < 10,          %if vector t is supplied together with a matrix y, a precise, time dependent guess will be constructed
    preciseguess=0;
else 
    preciseguess=1; 
end

if nargin < 9, 
    statev=0; 
    start=0;            %external function determines initial conditions because of no input for statev
    stateshock=0;       %no shock via state variables
elseif statev'*statev==0; 
    start=0;            %external function determines initial conditions because of input 'statev=0'
    stateshock=0;       %no shock via state variables
elseif sum((statev-1).^2)~=0
    stateshock=1;       %shock determined by change of state variables
    start=1;
else
    stateshock=0;       %shock determined by shock.m
    start=1;
end

if nargin < 8, M=50; end
if nargin < 7, y=ones(N,1); end
if nargin < 6, nu=0.1; end
if nargin < 5, n3=0; end
if nargin < 4, errorcode=1; end

%determination of guess for steady state values 
if preciseguess==0
    finalguess=y;
else
    finalguess=y(:,M);
end

ode=feval(funcODE,inf,finalguess);
if n3>0
    stat=feval(funcSTAT,inf,finalguess);
    res=[ode;stat]'*[ode;stat];
else
    res=ode'*ode;
    stat=[];
end


%Test, if funcODE and funcSTAT return real values
if ~(res<inf) | ~isreal(res)
    disp('ERROR: funcODE or funcSTAT return non real entries!');
    errorcode=2;
end

%Test, if funcODE and funcSTAT return vectors of right dimension
[test1 test2]=size(ode);
[test3 test4]=size(stat);
if ~((test1==n)&(test2==1)&(test3==n3)&(test4==1)) & n3>0
    disp('ERROR: funcODE or funcSTAT returns vector in wrong dimension!');
    errorcode=2;
elseif ~((test1==n)&(test2==1))
    disp('ERROR: funcODE returns vectors in wrong dimension!');
    errorcode=2;
end

if (sqrt(res) < 10^(-10)) & errorcode==0 %if guess for final steady state values is good enough...
    disp('Guess for steady state is sufficiently good.');
    disp(' ');
    if (stateshock==0) & (start==1)     %determine start via shock.m
        disp(['Calculation of initial steady state values...']);
        [start,fval,exitflag]=fsolve(@attach,finalguess,optimset('Display','off','MaxFunEvals',1000000,'MaxIter',100000),-1,funcODE,funcSTAT);
        if exitflag>0
            y=real(y);
            disp(['Convergence achieved.']);
            disp([' ']);
        else
            disp(['No convergence for initial steady state values!']);
            errorcode=3;
        end
    elseif stateshock==1 & start==1     %determine start via reduction in state variables
        start=finalguess.*statev;
    end                                 %else: external function determines initial conditions because of input 'statev=0'
elseif errorcode==0                     %guess for final steady state values is not good enough...
    if ~(preciseguess==1 & start==0)    %in the case of a precise guess with external initial conditions, nothing has to be calculated
        disp(['Calculation of final steady state values...']);
        [finalguess,fval,exitflag]=fsolve(@attach,finalguess,optimset('Display','off','MaxFunEvals',1000000,'MaxIter',100000),inf,funcODE,funcSTAT);
        if exitflag>0
            y=real(y);
            disp(['Convergence achieved.']);
            disp([' ']);
        else
            disp(['No convergence for final steady state values!']);
            errorcode=3;
        end
        if (stateshock==0) & (start==1) & (errorcode==0)    %determine start via shock.m
            disp(['Calculation of initial steady state values...']);
            [start,fval,exitflag]=fsolve(@attach,finalguess,optimset('Display','off','MaxFunEvals',1000000,'MaxIter',100000),-1,funcODE,funcSTAT);
            if exitflag>0
                y=real(y);
                disp(['Convergence achieved.']);
                disp([' ']);
            else
                disp(['No convergence for initial steady state values!']);
                errorcode=3;
            end
        elseif stateshock==1  & start==1     %determine start via reduction in state variables
            start=finalguess.*statev;
        end                                                 %else: determine initial conditions via external file
    end
end

guess=ones(M*N,1);
if (preciseguess==0) & (errorcode==0)    %As an initial guess the steady state values are used
    for i=1:M
        for ii=1:N
            guess((i-1)*N+ii)=finalguess(ii);
        end;
    end;
elseif errorcode==0
    tt=zeros(1,M);
    tau=[0:1/(M-1):1];
    for i=1:M
        if tau(i)~=1
            tt(i)=tau(i)/nu/(1-tau(i));
        else
            tt(i)=inf;
        end;
    end;
    y=releval(tt,t,y);
    for i=1:M
        guess((i-1)*N+1:i*N)=y(:,i);
    end
end

function attach=attach(x,t,funcODE,funcSTAT);
attach=[feval(funcODE,t,x);feval(funcSTAT,t,x)];
