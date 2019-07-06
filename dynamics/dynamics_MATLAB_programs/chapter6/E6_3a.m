% example 6.3(a)
% dynamics of a pendulum
% using an inline function

clear all; clc; close all
syms L m g t 

theta = sym('theta(t)');
omega = diff([0,0,theta],t);
alpha = diff(omega,t);
% diff(X,'t') or diff(X,sym('t')) 
% differentiates a symbolic expression 
% X with respect to t
% diff(X,'t',n) and diff(X,n,'t') 
% differentiates X n times 
% n is a positive integer

xC = L*cos(theta)/2;
yC = L*sin(theta)/2;

rC = [xC yC 0];
G = [0 -m*g 0];
IC = m*L^2/12;
IO = IC + m*(L/2)^2;
MO = cross(rC,G);
eq = -IO*alpha+MO;
eqz = eq(3);

eqI = subs(eqz,{L,m,g},{1,1,9.8});
eqI1 = subs(eqI,diff('theta(t)',t,2),'ddtheta');
eqI2 = subs(eqI1,diff('theta(t)',t),sym('x(2)'));
eqI3 = subs(eqI2,'theta(t)',sym('x(1)'));

% solve a second order ODE using MATLAB 
% ODE45 function
% write the second order equation as a 
% system of two first order equations, 
% by introducing x(2) = x(1)'
%  x(1)' = x(2)  ==>   x' = g(t,x)
%  x(2)' = f
%  
% define the vector x = [x(1); x(2)] 
% (column-vector). 

% first differential equation 
dx1 = sym('x(2)');
% second differential equation 
dx2 = solve(eqI3,'ddtheta');

% define right-hand side vector
dx1dt = char(dx1);
dx2dt = char(dx2);

% digits(3)
fprintf('equations of motion\n') 
fprintf('d[x(1)]/dt=%s\n',char(vpa(dx1dt))) 
fprintf('d[x(2)]/dt=%s\n',char(vpa(dx2dt))) 

eom=...
inline(sprintf('[%s; %s]',dx1dt,dx2dt),'t','x');
% inline(e) is an inline function object from the
% expression contained in the string 'e' 


t0 = 0;  % define initial time
tf = 5;  % define final time             
time = [0 tf];

x0 = [0; 0]; % define initial conditions

[t,xs] = ode45(eom, time, x0); 
% find t, xs, but don't show
% ode45  solves non-stiff differential equations, 
% medium order method.
% [ts,ys] = ode45(f,tspan,y0) integrates 
% the system of differential equations 
% y' = f(t,y) with initial conditions y0.     

x1 = xs(:,1); % extract x1 & x2 components from xs
x2 = xs(:,2);  

subplot(3,1,1),plot(t,x1,'r'),...
xlabel('t'),ylabel('\theta [rad]'),grid,...
subplot(3,1,2),plot(t,x2),...
xlabel('t'),ylabel('\omega [rad/s]'),grid,...
subplot(3,1,3),plot(x1,x2),...

xlabel('\theta [rad]'),ylabel('\omega [rad/s]'),grid

[ts,xs] = ode45(eom, 0:0.5:tf, x0);
format short
[ts, xs]

% end of program