% Simulation of stochastic control system
% To run the program pass the process noise covariance q=1,10,100, e.g. f6p1d4(10) 
function f6p1d4(q1)
    global u i gnoise q;
    q=q1;
    b=1; r=1; i=1; s=0; x1=[10 10]; % initial values 
    gnoise=randn(102,1); % white gaussian noise
    options = odeset('RelTol',1e-6,'AbsTol',[1e-6 1e-6]);
    for t=0:0.1:10
        [td s]=ode23(@sfunc,[t t+0.1],s);
        s=s(end);
        % Optimal Feedback gains & control input
        k=b*s/r; u(1)=-k*x1(end,1); u(2)=k; 
        [td x]=ode23(@xfunc,[t t+0.1],x1(end,:),options);
        x1(i,:)=x(end,:); t1(i)=td(end);
        i=i+1;
    end
    % Plotting
    plot(t1,x1(:,1),'k-',t1,x1(:,2),'k-.'); grid on; 
    legend('x(t)','X(t)'); title(strcat('q = ',num2str(q)));
    xlabel('Time (secs)'); axis([0 10 -5 12]);
return

% Function for Ricatti equation solution
function ds=sfunc(t,s)
    global q;
    a=0.05; b=1; r=1;
    ds(1)=2*a*s(1)-b^2*(s(1)^2)/r+q;
return

% Function for plant dynamics and lyapunov equation
function dx=xfunc(t,x)
    global u i gnoise;
    a=0.05; b=1; qp=5; g=1;
    dx(1)=a*x(1)+b*u(1)+g*qp*gnoise(i,1);
    dx(2)=2*(a-b*u(2))*x(2)+(g^2)*qp;
    dx=dx';
return