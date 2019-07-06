% example 7.3(b)
% Lagrange e.o.m.
clear all; clc; close all 

syms m l IC g real

t = sym('t','real');
q1 = sym('q1(t)');
q2 = sym('q2(t)');
q3 = sym('q3(t)');

rC = [q1 q2 0];
vC = diff(rC,t);
omega=[0 0 diff(q3,t)];

T=1/2*m*vC*vC.'+1/2*IC*omega*omega.';
T=simple(T);
fprintf('T = %s \n',char(T))

G = [0 -m*g 0]; 
rC_1 = deriv(rC, q1);
rC_2 = deriv(rC, q2);
rC_3 = deriv(rC, q3);
% generalized active forces
Q1 = rC_1*G.';
Q2 = rC_2*G.';
Q3 = rC_3*G.';

fprintf('Q1 = %s \n',char(Q1))
fprintf('Q2 = %s \n',char(Q2))
fprintf('Q3 = %s \n',char(Q3))
fprintf('\n')

fprintf('Lagrange e.o.m:\n')
%function [L] = Lagrange(T,Q,q,t)
Lagrange1 = Lagrange(T,Q1,q1,t);
Lagrange2 = Lagrange(T,Q2,q2,t);
Lagrange3 = Lagrange(T,Q3,q3,t);
fprintf(' %s =0 \n',char(Lagrange1))
fprintf(' %s =0 \n',char(Lagrange2))
fprintf(' %s =0 \n\n',char(Lagrange3))

% eom:
%  m*diff(q1(t), t, t) =0 
%  m*diff(q2(t), t, t) + g*m =0 
%  IC*diff(q3(t), t, t) =0 

% initial conditions at t=0
% q1(0)=0
% q2(0)=0
% q3(0)=q30
% Dq1(0)=v0*cos(alpha)
% Dq2(0)=v0*sin(alpha)
% Dq3(0)=IC*omega0

syms alpha q30 v0 omega0 real

q1f=dsolve...
('D2q1=0', 'q1(0)=0','Dq1(0)=v0*cos(alpha)');
q2f=dsolve...
('D2q2=-g','q2(0)=0','Dq2(0)=v0*sin(alpha)');
q3f=dsolve...
('D2q3=0', 'q3(0)=q30','Dq3(0)=omega0');

fprintf('q1(t) = %s \n',char(q1f))
fprintf('q2(t) = %s \n',char(q2f))
fprintf('q3(t) = %s \n',char(q3f))

% end of program

% q1(t) = t*v0*cos(alpha) 
% q2(t) = t*v0*sin(alpha) - (g*t^2)/2 
% q3(t) = q30 + omega0*t 
