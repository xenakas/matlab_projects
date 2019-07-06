% Example 4.2-1
% Optimal error covariance using Runge-Kutta integrator
function f4p2d1()
    global eps;
    eps_set=[0 0.3 0.7 1.0 1.5 2.0]; % Different values for epsilon
    for j=1:6
        eps=eps_set(j);
        [t p] = ode45(@errcov422,[0 3],[0 0 0]); % Runge-Kutta Integrator
        subplot(3,2,j); 
        plot(t,p(:,1),'k-',t,p(:,2),'k--',t,p(:,3),'k-.'); % Plot p
        axis([0 3 -0.1 0.5]); grid on; xlabel('t (secs)')
        title(strcat('\epsilon = ',num2str(eps)));
    end
return

% Differential Equations
function dp=errcov422(t,p)
    global eps
    dp = zeros(3,1);
    a1=-2; a2=-1; q1=1; q2=1; r1=1; r2=1;
    p1=p(1); p12=p(2); p2=p(3);
    dp(1)=-p1^2/r1-p12^2/r2+2*a1*p1+2*eps*p12+q1;
    dp(2)=-p1*p12/r1-p12*p2/r2+(a1+a2)*p12+eps*p2;
    dp(3)=-p12^2/r1-p2^2/r2+2*a2*p2+q2;
return