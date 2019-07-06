% example 4.7

clear all; clc; close all

syms m g R theta v0  N t 

theta = sym('theta(t)');
omega = diff(theta,t);
alpha = diff(omega,t);

% m at = - m g cos(theta)
% m an = m g sin(theta)+N
% =>
% m R theta'' = - m g cos(theta)
% m R (theta')^2 = m g sin(theta)+N
%
% eom:  
% R theta'' = - g cos(theta)  |*(theta')
% R (theta') theta'' = (theta') g cos(theta)

lhs= R*omega*alpha;
rhs=-omega*g*cos(theta);

ilhs=int(lhs);
irhs=int(rhs);

%int(R(theta')theta'')=int((theta') g cos(theta))+c

% t=0 => theta=0 and R*diff(theta,t)=v0 
cs = ilhs-irhs;
c = subs(cs,{diff(theta,t),theta},{v0/R,0});

eq = ilhs-irhs-c;

ql = {diff(theta,t,2),diff(theta,t),theta};   
qf = {'ddth', 'dth', 'th'};

eqs = subs(eq,diff(theta,t),'dth');
dtheta = solve(eqs,'dth');
dtheta = dtheta(1);

fprintf...
('equation of motion \n')
fprintf...
('d(theta)/dt = %s \n\n',char(dtheta))

N = m*R*(dtheta)^2- m*g*sin(theta);
N = simplify(N);
fprintf('N = %s \n\n',char(N))

syms theta0 real

slist={m, g, R, theta0, v0};
nlist={1, 10, 1, 0, 1};

omega0 = v0/R;
% m*R*alpha0 = - m*g*cos(theta0)
alpha0 = -g*cos(theta0)/R;
alphan = subs(alpha0,slist,nlist);
fprintf...
('alpha(0)=%s=%6.3f (rad/s^2)\n',char(alpha0),alphan)

% m R omega0^2 = m g sin(theta0)+N
N0 = m*R*omega0^2 - m*g*sin(theta0);
Nn = subs(N0,slist,nlist);
fprintf...
('N(0)=%s=%6.3f (N)\n',char(N0),Nn)

at = R*alpha0;
an = -R*omega0^2;
atn = subs(at,slist,nlist);
ann = subs(an,slist,nlist);
fprintf('at(0)=%6.3f (m/s^2)\n',atn)
fprintf('an(0)=%6.3f (m/s^2)\n',ann)
a = sqrt(atn^2+ann^2);
fprintf('a(0)=%6.3f (m/s^2)\n',a)

% end of program
%%

equation of motion 
d(theta)/dt = (v0^2 - 2*R*g*sin(theta(t)))^(1/2)/R 

N = (m*v0^2)/R - 3*g*m*sin(theta(t)) 

alpha(0)=-(g*cos(theta0))/R=-10.000 (rad/s^2)
N(0)=(m*v0^2)/R - g*m*sin(theta0)= 1.000 (N)
at(0)=-10.000 (m/s^2)
an(0)=-1.000 (m/s^2)
a(0)=10.050 (m/s^2)