% example 7.3(a)
% Hamilton e.o.m.
clear all; clc; close all 

syms m l IC g real

t = sym('t','real');
q1 = sym('q1(t)');
q2 = sym('q2(t)');
q3 = sym('q3(t)');
p1 = sym('p1(t)');
p2 = sym('p2(t)');
p3 = sym('p3(t)');

rC = [q1 q2 0];
vC = diff(rC,t);
omega=[0 0 diff(q3,t)];

T=m*vC*vC.'/2+IC*omega*omega.'/2;
T=simple(T);
fprintf('T = %s \n',char(T))
V = m*g*q2;

L = T - V;
p_1 = deriv(T, diff(q1,t));
p_2 = deriv(T, diff(q2,t));
p_3 = deriv(T, diff(q3,t));

fprintf('L = %s \n',char(L))
fprintf('generalized momenta: \n')
fprintf('p1 = %s \n',char(p_1))
fprintf('p2 = %s \n',char(p_2))
fprintf('p3 = %s \n\n',char(p_3))

% Hamiltonian of the system
H=p_1*diff(q1,t)+p_2*diff(q2,t)+p_3*diff(q3,t)-L;
fprintf('H = %s \n\n',char(H))

% Hamiltonian of the system expressed 
% in terms of generalized momenta
H=p_1*diff(q1,t)+p_2*diff(q2,t)+...
  p_3*diff(q3,t)-L;
H=subs(H,{diff(q1,t),diff(q2,t),diff(q3,t)},...
         {p1/m,p2/m,p3/IC});
H=expand(H);
          
fprintf('H = %s \n\n',char(H))

fprintf('Hamilton equations of motion:\n')
eq1 = diff(q1,t) - deriv(H, p1);
eq2 = diff(q2,t) - deriv(H, p2);
eq3 = diff(q3,t) - deriv(H, p3);
eq4 = diff(p1,t) + deriv(H, q1);
eq5 = diff(p2,t) + deriv(H, q2);
eq6 = diff(p3,t) + deriv(H, q3);

fprintf('%s = 0, (1)\n',char(eq1))
fprintf('%s = 0, (2)\n',char(eq2))
fprintf('%s = 0, (3)\n',char(eq3))
fprintf('%s = 0, (4)\n',char(eq4))
fprintf('%s = 0, (5)\n',char(eq5))
fprintf('%s = 0. (6)\n',char(eq6))
fprintf('\n')

% Hamilton equations of motion:
% diff(q1(t), t) - p1(t)/m = 0, (1)
% diff(q2(t), t) - p2(t)/m = 0, (2)
% diff(q3(t), t) - p3(t)/IC = 0, (3)
% diff(p1(t), t) = 0, (4)
% diff(p2(t), t) + g*m = 0, (5)
% diff(p3(t), t) = 0. (6)

% initial conditions t=0
% q1(0)=0
% q2(0)=0
% q3(0)=q30
% p1(0)=m*v0*cos(alpha)
% p2(0)=m*v0*sin(alpha)
% p3(0)=IC*omega0

syms q30 P1 P2 P3 real

sys=...
'Dq1=p1/m,Dq2=p2/m,Dq3=p3/IC,Dp1=0,Dp2=-m*g,Dp3=0';

ic=...
'q1(0)=0,q2(0)=0,q3(0)=q30,p1(0)=P1,p2(0)=P2,p3(0)=P3';

sol=dsolve(sys,ic);

q1t=sol.q1;
q2t=sol.q2;
q3t=sol.q3;
p1t=sol.p1;
p2t=sol.p2;
p3t=sol.p3;

syms alpha v0 omega0 real
plist ={P1,P2,P3};
plist0={m*v0*cos(alpha),m*v0*sin(alpha),IC*omega0};

q1f=subs(q1t,plist,plist0);
q2f=subs(q2t,plist,plist0);
q3f=subs(q3t,plist,plist0);
p1f=subs(p1t,plist,plist0);
p2f=subs(p2t,plist,plist0);
p3f=subs(p3t,plist,plist0);

fprintf('q1(t) = %s \n',char(q1f))
fprintf('q2(t) = %s \n',char(q2f))
fprintf('q3(t) = %s \n',char(q3f))
fprintf('p1(t) = %s \n',char(p1f))
fprintf('p2(t) = %s \n',char(p2f))
fprintf('p3(t) = %s \n',char(p3f))

% end of program

% q1(t) = t*v0*cos(alpha) 
% q2(t) = t*v0*sin(alpha) - (g*t^2)/2 
% q3(t) = q30 + omega0*t 
% p1(t) = m*v0*cos(alpha) 
% p2(t) = m*v0*sin(alpha) - g*m*t 
% p3(t) = IC*omega0 
