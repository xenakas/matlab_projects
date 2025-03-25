% Neoclassical Stochastic Growth Model
% Example based on Diaz-Gimenez (2001) - Linear Quadratic Approximations: An Introduction 
% (Chapter 2 in 'Computational Methods for the Study of Dynamic Economies', edited by Marimon & Scott)

% This model is also used by Aldrich, Fernandez-Villaverde, Gallant, & Rubio-Ramirez (2011) - "Tapping the supercomputer under your desk: Solving dynamic equilibrium models with graphics processors,"
% But they do use slightly different parameters to those used here.

Javier=0;   %If you set this parameter to 0 then the parameters will all be set to those used by Aldrich, Fernandez-Villaverde, Gallant, & Rubio-Ramirez (2011)

%% Set up
tauchenoptions.parallel=2; % Use GPU (is anyway the default option)
vfoptions.parallel=2; % Use GPU (is anyway the default option)

% The sizes of the grids
n_z=2^2;
n_k=2^12;

%Discounting rate
beta = 0.96;

%Give the parameter values (Params will be a 'structure' containing all the parameter values)
Params.alpha = 0.33;
Params.gamma=1; %gamma=1 is log-utility
Params.rho = 0.95;
Params.delta = 0.10;
Params.sigmasq_epsilon=0.09;

if Javier==0
    n_z=4;
    Params.beta=0.984;
    Params.gamma=2;
    Params.alpha=0.35;
    Params.delta=0.01;
    Params.rho=0.95;
    Params.sigma_epsilon=0.005;
    Params.sigmasq_epsilon=Params.sigma_epsilon^2;
    vfoptions.tolerance=(1-beta)*10^(-8); % Toolkit default is 10^(-9)
    vfoptions.howards=20; % Toolkit default is 80
end

% Params has been created as a structure. You can create the individual
% parameters from the structure by running the following command
CreateIndividualParams(Params)

%% Compute the steady state
K_ss=((alpha*beta)/(1-beta*(1-delta)))^(1/(1-alpha));
X_ss= delta*K_ss;
%These are not really needed; we just use them to determine the grid on
%capital. I mainly calculate them to stay true to original article.

%% Create grids (grids are defined as a column vectors)

q=3; % A parameter needed for the Tauchen Method
[z_grid, pi_z]=TauchenMethod(0,sigmasq_epsilon,rho,n_z,q,tauchenoptions); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q), transmatix is (z,zprime)

k_grid=linspace(0,20*K_ss,n_k)'; % Grids should always be declared as column vectors

%% Now, create the return function
DiscountFactorParamNames={'beta'};

ReturnFn=@(aprime_val, a_val, s_val, gamma, alpha, delta) StochasticNeoClassicalGrowthModel_ReturnFn(aprime_val, a_val, s_val, gamma, alpha, delta);
ReturnFnParamNames={'gamma', 'alpha', 'delta'}; %It is important that these are in same order as they appear in 'StochasticNeoClassicalGrowthModel_ReturnFn'

%% Solve
%Do the value function iteration. Returns both the value function itself,
%and the optimal policy function.
d_grid=0; %no d variable
n_d=0; %no d variable

tic;
V0=ones(n_k,n_z);
[V, Policy]=ValueFnIter_Case1(V0, n_d,n_k,n_z,d_grid,k_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
time=toc;

fprintf('Time to solve the value function iteration was %8.2f seconds. \n', time)


%% Draw a graph of the value function

%surf(V)

% Or to get 'nicer' x and y axes use
surf(k_grid*ones(1,n_z),ones(n_k,1)*z_grid',V)


