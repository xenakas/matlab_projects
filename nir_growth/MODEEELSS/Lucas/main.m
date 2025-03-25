% Version 3.1
% RELAXATION algorithm to solve infinite-horizon continuous time models.
%
% Description of procedure:
% Trimborn, Koch, Steger (2007) Multidimensional Transitional Dynamics: A
% Simple Numerical Procedure, forthcoming in Macroeconomic Dynamics
% 
% Copyright by Trimborn, Koch, Steger, 2008
%
% For further information contact Timo Trimborn, University of Hannover
% or visit http://www.relaxation.uni-siegen.de

clear all
disp(['Initialize Relaxation algorithm']);disp([' ']);

tic

globalpar                           % Initializes the global parameter

parini                              % Loads the Parameter Values

relaxsetting                        % Loads the settings for the Relaxation algorithm, i.e. dimensions, boundary conditions etc. 

% Converts settings to a form suitable for relax.m
[guess, start, errorcode]=initrelax(@funcODE, @funcSTAT, n, n1, n3, nu, y, M, statev);       

if errorcode==0            % Executes the relaxation algorithm if no error occured during initilization
    
    [t, x]=relax(@funcODE, @funcSTAT, @funcINI, @funcfinal, n, n1, n3, nu, guess, M, start, Endcond, maxit, tol, damp, dampfac);    
    
    %Normalization of specified variables
    for i=1:M
        x(normal,i)=x(normal,i)./x(normal,end);
    end;
    
    varex                               % Extracts the variables and stores them in the memory
 
%    If you want to calculate and display the eigenvalues at the steady
%    state, remove the comment of the two subsequent lines.
%    [EVa EVe Jac]=eigDAS(@funcODE, @funcSTAT, x(:,end));
%    disp(['Eigenvalues: ',num2str(EVa')]);disp([' ']);
    
end


%Calculation time
time=toc;
timesec=mod(time,60);
timemin=floor(time/60);
disp(['Calculation time: ',num2str(time),' seconds (',num2str(timemin),' min ',num2str(timesec),' sec)']);

%To get a first impression of the results, remove the comment of the
%subsequent line
plotrelax(t, x, n1, 100)
