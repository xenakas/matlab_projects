function varargout = relax2e(funcODE, funcSTAT, funcINI, funcfinal, n, n1, n3, nu, guess, M, start, Endcond, maxit, tol, damp, dampfac)
% Version 3.1
% RELAXATION algorithm to solve infinite-horizon continuous time models.
% This version includes an estimation of the global relative error along
% the mesh. The output solution is of second order, but obtained with the
% double number of mesh points than indicated by M.
%
% Description of procedure:
% Trimborn, Koch, Steger (2007) Multidimensional Transitional Dynamics: A
% Simple Numerical Procedure, forthcoming in Macroeconomic Dynamics
% 
% Copyright by Trimborn, Koch, Steger
%
% For further information contact Timo Trimborn, University of Hannover
% or visit http://www.relaxation.uni-siegen.de
%
% COMMAND:
% [t, x, error, errormean] = relax2e(@funcODE, @funcSTAT, @funcINI, @funcfinal, n, n1, n3, nu, y, M, start, Endcond, maxit, tol, damp, dampfac)
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
% error: (optional) relative error of the solution along the mesh (for each variable)
% errormean: (optional) relative mean error of the solution along the mesh

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


M=2*M-1;
y=ones(M*N,1);
for iii=1:2:M-1
    y((N*iii)-(N-1):(N*iii))=guess((N*(iii+1)/2)-(N-1):(N*(iii+1)/2));
    y((N*iii)+1:(N*iii)+N)=guess((N*(iii+1)/2)-(N-1):(N*(iii+1)/2));
end
y((N*M)-(N-1):(N*M))=guess((N*(M+1)/2)-(N-1):(N*(M+1)/2));

[t_d, x_d]=relax(funcODE, funcSTAT, funcINI, funcfinal, n, n1, n3, nu, y, M, start, Endcond, maxit, tol, damp, dampfac); 

M=1/2*(M+1);

guess(1:end)=x_d(:,1:2:2*M-1);

[t, x]=relax(funcODE, funcSTAT, funcINI, funcfinal, n, n1, n3, nu, guess, M, start, Endcond, maxit, tol, damp, dampfac);

%error of the solution of the dense mesh
x_error_d=1/3*(x_d(:,1:2:2*M-1)-x);
%calculatin of mean error along the mesh
for iii=1:M
    x_error_d_mean(iii)=sqrt(sum((x_error_d(:,iii)./x(:,iii)).^2));
end

varargout(1)={t};
varargout(2)={x_d(:,1:2:2*M-1)};
varargout(3)={x_error_d./x};
varargout(4)={x_error_d_mean};

