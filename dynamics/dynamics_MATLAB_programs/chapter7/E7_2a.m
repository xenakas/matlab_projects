% example 7.2.a
% Lagrange e.o.m.
clear all; clc; close all 

syms h alpha m1 m2 g real
      
t = sym('t','real');
q1 = sym('q1(t)');
q2 = sym('q2(t)');

r1=[q1 h/3 0];

a=h*cot(alpha)/3;
x2=q1+a-q2*cos(alpha);
y2=h-q2*sin(alpha);
r2=[x2 y2 0];

v1=diff(r1,t);
v2=diff(r2,t);

T1=1/2*m1*dotproduct(v1,v1);
T2=1/2*m2*dotproduct(v2,v2);
T = T1+T2;
fprintf('T1 = %s \n',char(T1))
fprintf('T2 = %s \n',char(T2))
fprintf('T = %s \n',char(T))
fprintf('\n')

% generalized active forces
G1 = [0 -m1*g 0]; 
G2 = [0 -m2*g 0];

r1_1 = deriv(r1, q1);
r1_2 = deriv(r1, q2);
r2_1 = deriv(r2, q1);
r2_2 = deriv(r2, q2);

% generalized active force Q1
Q1 = r1_1*G1.'+r2_1*G2.';
% generalized active force Q2
Q2 = r1_2*G1.'+r2_2*G2.';

fprintf('Q1 = %s \n',char(Q1))
fprintf('Q2 = %s \n',char(Q2))
fprintf('\n')

% another way of calculating Q1 and Q2
% potential energy V
V=m1*g*r1(2)+m2*g*r2(2);
Q_1=-deriv(V,q1);
Q_2=-deriv(V,q2);
fprintf('Q_1 = %d  \n',Q_1)
fprintf('Q_2 = %s  \n',char(Q_2))
fprintf('\n')

% first Lagrange's equation of motion 
% function [L] = Lagrange(T,Q,q,t)
Lagrange1 = Lagrange(T,Q1,q1,t);
% second Lagrange's equation of motion
Lagrange2 = Lagrange(T,Q2,q2,t);

Lagrange1 = simplify(Lagrange1);
Lagrange2 = simplify(Lagrange2);
fprintf('equations of motion:\n')
fprintf('%s =0  \n',char(Lagrange1))
fprintf('%s =0  \n',char(Lagrange2))
fprintf('\n')

slist = {m1, m2, g, alpha};
nlist = {10, 1, 9.81, pi/6};

Lagran1 = subs(Lagrange1, slist, nlist);
Lagran2 = subs(Lagrange2, slist, nlist);

digits(3)
fprintf('%s =0  \n',char(vpa(Lagran1)))
fprintf('%s =0  \n',char(vpa(Lagran2)))
fprintf('\n')

ql = {diff(q1,t,2), diff(q2,t,2)};    
qf = {'ddq1', 'ddq2'};

Lagra1 = subs(Lagran1, ql, qf);
Lagra2 = subs(Lagran2, ql, qf);

% solve e.o.m. for ddq1, ddq2
sol=solve(Lagra1,Lagra2,'ddq1,ddq2');
Lagr1 = sol.ddq1;
Lagr2 = sol.ddq2;

fprintf('a1 = %g (m/s^2)\n',double(Lagr1))
fprintf('a2 = %g (m/s^2)\n',double(Lagr2))

% end of program


% equations of motion:
% m1*diff(q1(t), t, t) - (m2*(2*cos(alpha)*diff(q2(t), t, t) - 2*diff(q1(t), t, t)))/2 =0  
% m2*diff(q2(t), t, t) - g*m2*sin(alpha) - m2*cos(alpha)*diff(q1(t), t, t) =0  
%
% 11.0*diff(q1(t), t, t) - 0.866*diff(q2(t), t, t) =0  
% diff(q2(t), t, t) - 0.866*diff(q1(t), t, t) - 4.91 =0  
% 
% a1 = 0.414 (m/s^2)
% a2 = 5.26 (m/s^2)
