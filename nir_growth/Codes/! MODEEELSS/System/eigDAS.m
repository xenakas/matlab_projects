function varargout=eigDAS(funcODE, funcSTAT, y, time);
%Form: [EVa EVe Jac]=eigDAS(@funcODE, @funcSTAT, n_d, n_s, y);
%Calculates the eigenvalues of the linerized differential algebraic system 
%funcODE and funcSTAT at the point y numerically. funcODE and funcSTAT denote the 
%function handles of the differential equations and the algebraic equations, respectively.
%They are used according to funcODE(time,y), and funcSTAT(time,y), where time is the 
%point of time at which the output is calculated. Input of time is
%optional, default is time=inf.
%Note, that the system has to be of differential index one.
%First output: vector of eigenvalues
%Second output: matrix with column vectors of the corresponding
%eigenvectors
%Third output: Jacobian matrix
%Example: Call eigDAS(@funcODE, @funcSTAT, x(:,end));
%after you ran main.m of the relaxation algorithm to display the
%eigenvalues at the stationary point
%
%Copyright: Timo Trimborn, University of Hannover, 2008

if nargin < 4
    time=inf;
end 

n=length(funcODE(time,y));
n3=length(funcSTAT(time,y));
N=n+n3;

if length(y)~=n+n3
    disp(['ERROR: dim(y) does not equal dim(funcODE) + dim(funcSTAT)!']);
    varargout(1)=[];
    varargout(2)=[];
    varargout(3)=[];
else
    fac=[];
    [Jac1,fac]=numjac(funcODE,time,y,feval(funcODE,time,y),sqrt(eps)*ones(N,1),fac,0,[],[]);
    if n3>0
        [Jac2,fac]=numjac(funcSTAT,time,y,feval(funcSTAT,time,y),sqrt(eps)*ones(N,1),fac,0,[],[]);
        kernel=null(Jac2);
        CondNumb=cond(kernel(1:n,:));
        if CondNumb>1/sqrt(eps)
            disp(['ERROR: The linearized Differential Algebraic System is of higher Differential Index!']);
            disp(['Eigenvalues and Eigenvectors cannot be computed!']);
        end
        ReducedDS=inv(kernel(1:n,:))*Jac1*kernel;
        A=ReducedDS;
    else
        A=Jac1;
        ReducedDS=Jac1;
    end %if n3>0
    [EVe EVa]=eig(ReducedDS);
    if n3>0
        EVe=kernel*EVe;
    end
    varargout(1)={diag(EVa)};
    varargout(2)={EVe};
    varargout(3)={A};
end     %if length...



