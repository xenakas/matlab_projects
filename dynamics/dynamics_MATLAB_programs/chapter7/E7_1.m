% example 7.1
clear all; clc; close all 

syms L m M k g t real
q1 = sym('q1(t)');
q2 = sym('q2(t)');

rA= [q1,0,0];
rB= [q1+L*sin(q2),L*cos(q2),0];

vA=diff(rA,t);
vB=diff(rB,t);

% kinetic energy of 1
T1 = m*vA*vA.'/2;
% .' array transpose 
% A.' is the array transpose of A
% kinetic energy of 2
T2 = M*vB*vB.'/2;
T2 = simplify(T2);
fprintf('T1 = %s \n',char(T1))
fprintf('T2 = %s \n',char(T2))
% total kinetic energy
T = expand(T1 + T2);
fprintf('T = \n'); pretty(simple(T)); 
fprintf('\n')

%deriv(f,g(t)) 
%differentiates f with respect to g(t)

% dT/d(dq)
Tdq1 = deriv(T, diff(q1,t));
Tdq2 = deriv(T, diff(q2,t));
fprintf('dT/d(dq1) = \n'); pretty(simple(Tdq1));
fprintf('\n')
fprintf('dT/d(dq2) = \n'); pretty(simple(Tdq2)); 
fprintf('\n')

% d(dT/d(dq))/dt
Tt1 = diff(Tdq1, t);
Tt2 = diff(Tdq2, t);
fprintf('d dT/d(dq1)/dt = \n'); 
pretty(simple(Tt1)); 
fprintf('\n')
fprintf('d dT/d(dq2)/dt = \n'); 
pretty(simple(Tt2)); 
fprintf('\n')

% dT/dq
Tq1 = deriv(T, q1);
Tq2 = deriv(T, q2);
fprintf('dT/dq1 = %g \n',Tq1)
fprintf('\n')
fprintf('dT/dq1 = %s \n',char(Tq2)) 
fprintf('\n')

% left hand side of Lagrange's eom
LHS1 = Tt1 - Tq1;
LHS2 = Tt2 - Tq2;

% generalized active forces
G1 = [0 m*g 0]; 
G2 = [0 M*g 0];
Fe1=[-k*q1 0 0];

% partial derivatives
rA_1 = deriv(rA, q1); rB_1 = deriv(rB, q1);
rA_2 = deriv(rA, q2); rB_2 = deriv(rB, q2);

% generalized active force Q1
Q1 = rA_1*(G1+Fe1).'+rB_1*G2.';
% generalized active force Q2
Q2 = rA_2*(G1+Fe1).'+rB_2*G2.';
fprintf('Q1 = \n'); pretty(simple(Q1)); 
fprintf('\n')
fprintf('Q2 = \n'); pretty(simple(Q2)); 
fprintf('\n')

% first Lagrange's equation of motion
Lagrange1 = LHS1-Q1;
% second Lagrange's equation of motion
Lagrange2 = LHS2-Q2;

data = {m, M, k, L, g};
datn = {0.4, 0.2, 0.75, 0.3, 9.81};

Lagran1 = subs(Lagrange1, data, datn);
Lagran2 = subs(Lagrange2, data, datn);

ql = {diff(q1,t,2), diff(q2,t,2),...
    diff(q1,t), diff(q2,t), q1, q2};    
qf = ...
{'ddq1', 'ddq2', 'x(2)', 'x(4)', 'x(1)', 'x(3)'};

% ql                    qf
%----------------------------
% diff('q1(t)',t,2) -> 'ddq1'
% diff('q2(t)',t,2) -> 'ddq2'
%   diff('q1(t)',t) -> 'x(2)'
%   diff('q2(t)',t) -> 'x(4)'
%           'q1(t)' -> 'x(1)'
%           'q2(t)' -> 'x(3)' 


Lagra1 = subs(Lagran1, ql, qf);
Lagra2 = subs(Lagran2, ql, qf);

% solve e.o.m. for ddq1, ddq2
sol = solve(Lagra1,Lagra2,'ddq1, ddq2');
Lagr1 = sol.ddq1;
Lagr2 = sol.ddq2;

% system of ODE 
dx2dt = char(Lagr1); 
dx4dt = char(Lagr2);

fid = fopen('eomE7_1.m','w+'); 
fprintf(fid,'function dx = eomE7_1(t,x)\n');
fprintf(fid,'dx = zeros(4,1);\n');
fprintf(fid,'dx(1) = x(2);\n');
fprintf(fid,'dx(2) = ');
fprintf(fid,dx2dt);
fprintf(fid,';\n');
fprintf(fid,'dx(3) = x(4);\n');
fprintf(fid,'dx(4) = ');
fprintf(fid,dx4dt);
fprintf(fid,';');
fclose(fid); cd(pwd);

t0 = 0;  tf = 5; time = [0 tf];

x0 = [0.05 0.1 pi/4 0]; % initial conditions

[t,xs] = ode45(@eomE7_1, time, x0);

x1 = xs(:,1); 
x2 = xs(:,2);
x3 = xs(:,3); 
x4 = xs(:,4);    

subplot(2,1,1),plot(t,x1,'r'),...
xlabel('t (s)'),ylabel('q1 (m)'),grid,...
subplot(2,1,2),plot(t,x3*180/pi,'b'),...
xlabel('t (s)'),ylabel('q2 (deg)'),grid

[ts,xs] = ode45(@eomE7_1,0:1:5,x0);
fprintf('Results \n\n')
fprintf...
('    t(s)      q1(m)   q2(rad)\n') 
[ts, xs(:,1), xs(:,3)]

% end of program