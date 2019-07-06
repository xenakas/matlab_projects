% Example 3.2-3
function f3p2d6()
    % Runge Kutta Integrator for Continuous time Kalman Filter
    options = odeset('RelTol',1e-4,'AbsTol',[1e-8 1e-8 1e-8]);
    % Runge-Kutta Integrator
    [t p] = ode45(@errcov,[0 10],[0 0 0],options);
    % Plot the error covariances
    plot(t,p(:,1),t,p(:,2),t,p(:,3));
    legend('p_1','p_2','p_3'); grid on; xlabel('t (secs.)');
return

function dp=errcov(t,p)
    dp = zeros(3,1);
    q=1; r=1; wns=0.64; a=0.16; % Constants
    dp(1)=2*p(2)-(1/r)*p(2)^2; % Differential Equations
    dp(2)=p(3)-wns*p(1)-2*a*p(2)-p(2)*p(3)/r;
    dp(3)=-2*wns*p(2)-4*a*p(3)-(1/r)*p(3).^2+q;
return