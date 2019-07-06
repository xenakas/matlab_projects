% example 7.2(b)
% Hamilton e.o.m.
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

% kinetic energy
T1=1/2*m1*v1*v1.';
T2=1/2*m2*v2*v2.';
T = T1+T2;
% potential energy V
V=m1*g*r1(2)+m2*g*r2(2);
% Lagrangian L
L = T-V;

p1L = deriv(L, diff(q1,t));
p2L = deriv(L, diff(q2,t));
fprintf('p1 = %s \n',char(simple(p1L)))
fprintf('p2 = %s \n',char(simple(p2L)))
fprintf('\n')

p1 = sym('p1(t)');
p2 = sym('p2(t)');

p1dq = subs...
(p1L,{diff(q1,t),diff(q2,t)},{'dq1','dq2'});
p2dq = subs...
(p2L,{diff(q1,t),diff(q2,t)},{'dq1','dq2'});

solp=solve(p1dq-p1,p2dq-p2,'dq1','dq2');

dq1p=simplify(solp.dq1);
dq2p=simplify(solp.dq2);

fprintf('dq1/dt = %s \n',char(dq1p))
fprintf('dq2/dt = %s \n',char(dq2p))
fprintf('\n')

% Hamiltonian H 
% H=sum(p_i*diff(q_i)) - L
Hdq=p1L*diff(q1,t)+p2L*diff(q2,t)-L;
Hdq=simplify(Hdq);        
fprintf('H = %s \n',char(simple(Hdq)))
fprintf('\n')

Hdp=simplify(subs...
(Hdq,{diff(q1,t),diff(q2,t)},{dq1p,dq2p}));
fprintf('H = %s \n',char(simple(Hdp)))
fprintf('\n')

% Hamiltonian equations 
Hdq1 =  deriv(Hdp, p1);
Hdq2 =  deriv(Hdp, p2);
Hdp1 = -deriv(Hdp, q1);
Hdp2 = -deriv(Hdp, q2);

fprintf('Hamilton eom \n\n')
fprintf('dq1/dt = %s \n',char(Hdq1))
fprintf('dq2/dt = %s \n',char(Hdq2))
fprintf('\n')
fprintf('dp1/dt = %d \n',Hdp1)
fprintf('dp2/dt = %s \n',char(Hdp2))
fprintf('\n')

dp1dt = diff(p1L,t);
dp2dt = diff(p2L,t);

eq1=dp1dt-Hdp1;
eq2=dp2dt-Hdp2;

fprintf('from Hamilton eom => \n')
fprintf('%s = 0 \n',char(eq1))
fprintf('%s = 0 \n',char(eq2))
fprintf('\n')

slist = {m1, m2, g, alpha};
nlist = {10, 1, 9.81, pi/6};

e1 = subs(eq1, slist, nlist);
e2 = subs(eq2, slist, nlist);

digits(3)
fprintf('%s =0  \n',char(vpa(e1)))
fprintf('%s =0  \n',char(vpa(e2)))
fprintf('\n')

ql = {diff(q1,t,2), diff(q2,t,2)};    
qf = {'ddq1', 'ddq2'};

h1 = subs(e1, ql, qf);
h2 = subs(e2, ql, qf);

% solve e.o.m. for ddq1, ddq2
sol=solve(h1,h2,'ddq1,ddq2');
a1 = sol.ddq1;
a2 = sol.ddq2;

fprintf('a1 = %g (m/s^2)\n',double(a1))
fprintf('a2 = %g (m/s^2)\n',double(a2))

% end of program

% 11.0*diff(q1(t), t, t) - 0.866*diff(q2(t), t, t) =0  
% 1.0*diff(q2(t), t, t) - 0.866*diff(q1(t), t, t) - 4.91 =0  
% 
% a1 = 0.414 (m/s^2)
% a2 = 5.26 (m/s^2)
