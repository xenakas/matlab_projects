%% Waveform Relaxation implementation for a Real Business Cycle model
% backward integration implementation of Brunner and Strulik (2002)

% clean up
clear all
warning('off', 'MATLAB:nearlySingularMatrix')
% close all
home

% fzero options
options=optimset('Display','off','MaxFunEvals',1000000,'MaxIter',100000,'TolFun',10^-15,'TolX',10^-15);
conv_crit = 10^8;                   % set initial convergence criteria           
conv_crit_old=conv_crit*2;
% start counting computation time
tic

%% continuous-time growth model, infinite horizon problem
% variables
% c(t)      : consumption (continuous control on R+)
% k(t)      : capital stock (continuous state on R+)
% y(t)      : aggregate output
% r(t)      : real interest rate on physical capital (net)
% w(t)      : real wage rate

%% select model parameters
global alpha delta rho theta css kss lambda gamma determ policy_l policy_u gamma2 alpha2 c_tilde_spline c_spline   spline_on direction kinit  c_tilde_spline 

alpha  = 0.5;     % output elasticity of capital
delta  = 0.05;    % rate of physical capital depreciation
theta  = 2.5;     % coefficient of relative risk aversion
lambda = 0.2;    % arrival rate disasters
gamma  = 0.1;    % size of disasters
% select rho such that the saving rate is constant (if applicable)
rho = ((1-gamma)^(1-alpha*theta)-1)*lambda-(1-alpha*theta)*delta; 
% rho    = 0.017838019216068;     % subjective rate of time preference

gamma2=gamma;
alpha2=alpha;

%% select options
% optional settings
spline_on = 1;                      % use spline interpolation
plots_on = 1;                       % displays various plots
figps_on = 0;                       % output of plots as postscript
i_max = 30;                         % maximium iterations

%% compute non-stochastic steady-state values
kss = (alpha/(rho+delta))^(1/(1-alpha)); 
css = kss^alpha - delta*kss;
rss = alpha*kss^(alpha-1)-delta;

%% select algorithm settings

% algorithm settings
tolerance_ODE = 10^-14;             % smaller values decrease approximation error (we use 10^-14 in the paper)
tol=10^-8;                          % convergence tolerance
% smooth=0.5;                         % smoothing parameter in csaps which determines the relative weight you would like to place on the contradictory demands of having f be smooth vs having f be close to the data
                                    % smooth = 0, f is the least-squares straight line fit to the data, while, at the other extreme, i.e., for smooth = 1, f is the variational, or `natural' cubic spline interpolant
                                    % smooth moves from 0 to 1, the smoothing spline changes from one extreme to the other

%% solve for the optimal policy function (lower arm)
direction=0;    % direction = 0 (no transformation)
determ=1;       % selects the deterministic system of ODEs

disp(['solve for the policy function (lower arm)']);
k_min=0;
kinit=kss-1e-12;
cinit=css;

[t,x,step] = ode45mod(@funcODE,[kinit;cinit],k_min,1,tolerance_ODE,0);
policy_l = x';

kss_old=kss;
css_old=css;

%% solve for the optimal policy function (upper arm)
disp(['solve for the policy function (upper arm)']);
direction=1;   % direction = 1 (transformation)

k_max = inf;
kinit=kss+1e-12;
cinit=css;

[t,x,step] = ode45mod(@funcODE,[1/kinit;1/cinit],1/k_max,1,tolerance_ODE,0);
policy_u = x(length(x'):-1:1,:)';

%% clean up and interpolate policy function and slope

% interpolation of the policy functions lower arm, upper arm, complete
if spline_on
    c_spline   = csapi([0 policy_l(1,:) 1./policy_u(1,:)],[0 policy_l(2,:) 1./policy_u(2,:)]);
else
    c_spline   = interp1([0 policy_l(1,:) 1./policy_u(1,:)],[0 policy_l(2,:) 1./policy_u(2,:)],'linear','pp');
end

c_tilde_l = ppval(c_spline,(1-gamma2)*policy_l(1,:))./ppval(c_spline,policy_l(1,:));
c_tilde_u = ppval(c_spline,(1-gamma2)*1./policy_u(1,:))./ppval(c_spline,1./policy_u(1,:));

% computing jump terms from the complete interpolation
if spline_on
    c_tilde_spline = csapi([policy_l(1,:) 1./policy_u(1,:)],[c_tilde_l c_tilde_u]);
else
    c_tilde_spline = interp1([policy_l(1,:) 1./policy_u(1,:)],[c_tilde_l c_tilde_u],'linear','pp');
end

% save policy function of deterministic case
policy_determ=[policy_l 1./policy_u];
kss_determ=kss;
css_determ=css;

%% Waveform Relaxation

determ=2;       % selects the conditional deterministic system of ODEs
index=0;
while conv_crit>tol  & index<i_max & conv_crit<conv_crit_old
    index=index+1;
    conv_crit_old=conv_crit;
    c_spline_last_iteration = c_spline;
    c_tilde_spline_last_iteration = c_tilde_spline;
    policy_last_iteration_l=policy_l;
    policy_last_iteration_u=policy_u;
    kss_last_iteration=kss;
    css_last_iteration=css;
    
    %--------------------------------------------------------------------------
    % Step 3,4   compute jump terms, find new steady state
    %--------------------------------------------------------------------------
    [kss,fval,exitflag]=fzero(@findss,kss,options);
    css=kss^alpha-delta*kss;
    rss = alpha*kss^(alpha-1)-delta;
    
    %--------------------------------------------------------------------------
    % Step 5   solve with jump sensitivity
    %--------------------------------------------------------------------------
    disp(['Iteration ',num2str(index)]);
    
    direction=0;    % direction = 0 (no transformation)
    disp(['solve for the policy function (lower arm)']);
    kinit=kss-1e-12;
    cinit=css;

    [t,x,step] = ode45mod(@funcODE,[kinit;cinit],k_min,1,tolerance_ODE,0);
    policy_l = x(x(:,1)'>0,:)';
    
    direction=1;   % direction = 1 (transformation)
    disp(['solve for the policy function (upper arm)']);
    kinit=kss+1e-12;
    cinit=css;
   
    [t,x,step] = ode45mod(@funcODE,[1/kinit;1/cinit],1/k_max,1,tolerance_ODE,0);
    policy_u = x(length(x'):-1:1,:)';


    % interpolate the policy function of the last iteration along the mesh
    policy_compare_last_iteration_l = ppval(c_spline_last_iteration, policy_l(1,:));
    policy_compare_last_iteration_u = 1./ppval(c_spline_last_iteration, 1./policy_u(1,:));
    
    if spline_on
        c_spline   = csapi([0 policy_l(1,:) 1./policy_u(1,:)],[0 policy_l(2,:) 1./policy_u(2,:)]);
    else
        c_spline   = interp1([0 policy_l(1,:) 1./policy_u(1,:)],[0 policy_l(2,:) 1./policy_u(2,:)],'linear','pp');
    end

    c_tilde_l = ppval(c_spline,(1-gamma2)*policy_l(1,:))./ppval(c_spline,policy_l(1,:));
    c_tilde_u = ppval(c_spline,(1-gamma2)*1./policy_u(1,:))./ppval(c_spline,1./policy_u(1,:));

    % computing jump terms from the complete interpolation
    if spline_on
        c_tilde_spline = csapi([policy_l(1,:) 1./policy_u(1,:)],[c_tilde_l c_tilde_u]);
    else
        c_tilde_spline = interp1([policy_l(1,:) 1./policy_u(1,:)],[c_tilde_l c_tilde_u],'linear','pp');
    end

    %----------------------
    %convergence criterium: deviation of the policy function between two
    %iterations (absolute and relative deviation)
    %----------------------



    % calculate the analytic policy function and convergence criterium, if applicable
 
    if (alpha ~= theta) 
        policy_true_l = (theta-1)/theta*abs(policy_l(1,:)).^alpha;
        policy_true_u = 1./((theta-1)/theta*abs(1./policy_u(1,:)).^alpha);
    else
        policy_true_l = (rho-(exp((1-theta)*log(1-gamma))-1)*lambda+(1-theta)*delta)/theta*abs(policy_l(1,:));
        policy_true_u = 1./((rho-(exp((1-theta)*log(1-gamma))-1)*lambda+(1-theta)*delta)/theta*abs(1./policy_u(1,:)));
    end
    

    state_eval=[0.5*kss:0.01:kss*2];

    conv_crit=norm(ppval(c_spline_last_iteration,state_eval)-ppval(c_spline,state_eval));
    disp(['convergence criterium (absolute deviation): ',num2str(conv_crit)]);

      
    % illustrates the results
    if plots_on
        % compare deterministic policy function with the stochastic one
        figure(101);
        plot_policy
        % compare jump terms
        figure(102);
        plot_ctilde
        % plot true errors, if applicable
        if alpha==theta || rho==((1-gamma)^(1-alpha*theta)-1)*lambda-(1-alpha*theta)*delta
            figure(104);
            plot_errors
            figure(105);
            plot_errors_rel
        end
        figure(106);
        plot_errors_last
         figure(107);
        plot_errors_last_rel
        pause(1);
    end

end %for


    
%Calculation time
time=toc;
timesec=mod(time,60);
timemin=floor(time/60);
disp(['Calculation time: ',num2str(time),' seconds (',num2str(timemin),' min ',num2str(timesec),' sec)']);
 

