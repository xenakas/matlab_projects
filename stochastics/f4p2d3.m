% Code to find the predicted and actual error covariances using 
% Runge-Kutta Integrator
function f4p2d3()
    global eps;
    eps_set=[0.3 0.7 1.0 1.5]; % Different values for epsilon
    for j=1:4
        eps=eps_set(j);
        [t s] = ode45(@s422,[0 3],[0 0 0 0 0 0 0 0 0 0 0 0 0]); % Runge-Kutta Integrator
        subplot(2,2,j); 
        plot(t,s(:,11),'k-',t,s(:,12),'k--',t,s(:,13),'k--',t,s(:,1),'k-.');
        axis([0 3 -0.1 0.5]); grid on; xlabel('t (secs)')
        title(strcat('\epsilon = ',num2str(eps)));
    end
return

% Function for actual error covariance
function ds=s422(t,s)

    ds = zeros(12,1);
    a1=-2; a2=-1; q1=1; q2=1; r1=1; r2=1;
    global eps

    % suboptimal (predicted) error covariance
    p1=s(1); p12=s(2); p2=s(3);
    ds(1)=-p1^2/r1-p12^2/r2+2*a1*p1+2*eps*p12+q1;
    ds(2)=-p1*p12/r1-p12*p2/r2+(a1+a2)*p12+eps*p2; 
    ds(3)=-p12^2/r1-p2^2/r2+2*a2*p2+q2; 

    % Suboptimal Gains
    K1=p1/r1;
    K2=p2/r2;

    % Plant state correlations
    u1=s(4); u12=s(5); u2=s(6);
    ds(4)=2*a1*u1+2*eps*u12+q1;
    ds(5)=(a1+a2)*u12+eps*u2;
    ds(6)=2*a2*u2+q2;

    % Actual error/plant state crosscorrelations 
    v1=s(7); v12=s(8); v21=s(9); v2=s(10);
    ds(7)=(2*a1-K1)*v12+eps*u12;
    ds(8)=(a1+a2-K1)*v12+eps*u2;
    ds(9)=(a1+a2-K2)*v21+eps*v2;
    ds(10)=(2*a2-K2)*v2+q2;

    % Actual error correlations
    s1=s(11); s12=s(12); s2=s(13);
    ds(11)=2*(a1-K1)*s1+r1*K1^2+2*eps*v12+q1;
    ds(12)=(a1+a2-K1-K2)*s12+eps*v2;
    ds(13)=2*(a2-K2)*s2+r2*K2^2+q2;
    
return