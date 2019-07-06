% example 6.5(a)
% RT chain
% Newton-Euler equations of motion

clear all; clc; close all;

L = 1;    % m
m = 1;    % kg
g = 9.81; % m/s^2
IA = 1;   % kg m^2

syms t

xC = L*cos(sym('theta(t)'))/2;
yC = L*sin(sym('theta(t)'))/2;
rC = [xC yC 0];

xA = sym('r(t)')*cos(sym('theta(t)'));
yA = sym('r(t)')*sin(sym('theta(t)'));
rA = [xA yA 0];

vA = diff(rA, t); 
aA = diff(vA, t);

omega = [0 0 diff('theta(t)', t)];
alpha = diff(omega, t);

syms f21 
F21x = -f21*sin(sym('theta(t)'));
F21y =  f21*cos(sym('theta(t)'));
F21 = [F21x F21y 0];

syms p
xP = p*cos(sym('theta(t)'));
yP = p*sin(sym('theta(t)'));
rP = [xP yP 0];

G1 = [0 -m*g 0];
G2 = [0 -m*g 0];

IO = m*L^2/3;
% - IO alpha + rC x G1 + rP x F21 = 0 (1)
e1 = -IO*alpha + cross(rC, G1) + cross(rP, F21);
e1z = e1(3);

% m2 aA  = - F21 + G2 => "
e2 = -m*aA - F21 + G2;
% (x):   m2 aAx + F21x = 0 (2)
e2x = e2(1);
% (y):   m2 aAy + F21y - G2 = 0 (3)
e2y = e2(2);

% -IA alpha + (rP-rA) x (-F21) = 0 (4)
e2A = -IA*alpha+cross(rP-rA,-F21);
e2Az = e2A(3);

sol = solve(e2y, e2Az, f21, p);
f21s = sol.f21;
ps = sol.p;

EqI  = subs(e1z,{'f21','p'},{f21s,ps});

EqII = subs(e2x,{'f21','p'},{f21s,ps});

% list for the symbolical variables 
slist={diff('theta(t)',t,2),diff('r(t)',t,2),...
diff('theta(t)',t),diff('r',t),'theta(t)','r(t)'};
nlist = ...
{'ddq1', 'ddq2', 'x(2)', 'x(4)', 'x(1)','x(3)'};
% diff('theta',t,2) will be replaced by 'ddq1'
% diff('r(t)',t,2) will be replaced by 'ddq2'
% diff('theta',t) will be replaced by 'x(2)'
% diff('r(t)',t) will be replaced by 'x(4)'
% 'theta' will be replaced by 'x(1)'
% 'r(t)' will be replaced by 'x(3)'

eq1 = subs(EqI,slist,nlist);
eq2 = subs(EqII,slist,nlist);

sole = solve(eq1,eq2,'ddq1, ddq2');
ddq1s = sole.ddq1;
ddq2s = sole.ddq2;

dx2dt = char(ddq1s);
dx4dt = char(ddq2s);

fid = fopen('eomE6_5a.m','w+'); 
fprintf(fid,'function dx = eomE6_5a(t,x)\n');
fprintf(fid,'dx = zeros(4,1);\n');
fprintf(fid,'dx(1) = x(2);\n');
fprintf(fid,'dx(2) = ');
fprintf(fid,dx2dt);
fprintf(fid,';\n');
fprintf(fid,'dx(3) = x(4);\n');
fprintf(fid,'dx(4) = ');
fprintf(fid,dx4dt);
fprintf(fid,';');
fclose(fid);
cd(pwd);

x0 = [pi/6 0 0.5 0]; % initial conditions

option = odeset...
('RelTol',1e-3,'MaxStep',1e-3,'Events',@eventE6_5);
t0 = 0;  tf = .5; time = [0 tf];

[t, xs, te, ye] = ode45(@eomE6_5a, time, x0, option); 

x1 = xs(:,1); 
x2 = xs(:,2);
x3 = xs(:,3); 
x4 = xs(:,4);  

subplot(2,1,1),plot(t,x1*180/pi,'r'),...
xlabel('t (s)'),ylabel('\theta (deg)'),grid,...
subplot(2,1,2),plot(t,x3,'b'),...
xlabel('t (s)'),ylabel('r (m)'),grid

[ts,xs] = ode45(@eomE6_5a,0:0.1:1,x0); 
[ts,xs]

% end of program