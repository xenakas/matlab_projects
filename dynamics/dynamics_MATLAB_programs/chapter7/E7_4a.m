% example 7.4
% Lagrange e.o.m.
clear all; clc; close all 

syms t L1 L2 L m1 m2 m g 

q1 = sym('q1(t)');
q2 = sym('q2(t)');

c1 = cos(q1); s1 = sin(q1);
c2 = cos(q2); s2 = sin(q2);

L1 = L; L2 = L;

xA = L1*s1; 
yA = L1*c1; 
rA = [xA yA 0];
rC1 = rA/2; 
vC1 = diff(rC1,t);
xB = xA + L2*s2; 
yB = yA + L2*c2; 
rB = [xB yB 0];
rC2 = (rA + rB)/2; 
vC2 = diff(rC2,t);
omega1 = [0 0 -diff(q1,t)]; 
omega2 = [0 0 -diff(q2,t)]; 

m1 = m; m2 = m;
IO = m1*L1^2/3; IC2 = m2*L2^2/12;

% kinetic energy of the link 1
T1 = IO*omega1*omega1.'/2;
fprintf('T1 = \n')
pretty(T1); fprintf('\n')

% kinetic energy of the link 2
T2 = m2*vC2*vC2.'/2 + IC2*omega2*omega2.'/2;
T2 = simplify(T2);
fprintf('T2 = \n')
pretty(T2); fprintf('\n')

% total kinetic energy
T = simplify(T1 + T2);
fprintf('T = \n') 
pretty(T); fprintf('\n')

%deriv(f,g(t)) differentiates f with respect to g(t)

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
fprintf('dT/dq1 = \n'); pretty(Tq1); 
fprintf('\n')
fprintf('dT/dq2 = \n'); pretty(Tq2); 
fprintf('\n')

% left hand side of Lagrange's eom
LHS1 = Tt1 - Tq1;
LHS2 = Tt2 - Tq2;

% generalized active forces
G1 = [0 m1*g 0]; 
G2 = [0 m2*g 0];

% partial derivatives
rC1_1 = deriv(rC1, q1); 
rC2_1 = deriv(rC2, q1);
rC1_2 = deriv(rC1, q2); 
rC2_2 = deriv(rC2, q2);

% generalized active force Q1
Q1 = rC1_1*G1.'+ rC2_1*G2.';
% generalized active force Q2
Q2 = rC1_2*G1.'+ rC2_2*G2.';

fprintf('Q1 = \n'); pretty(simple(Q1)); 
fprintf('\n')
fprintf('Q2 = \n'); pretty(simple(Q2)); 
fprintf('\n')


% another way of calculating Q1 and Q2
% potential energy V
V=-m1*g*rC1(2)-m2*g*rC2(2);
Q_1=-deriv(V,q1);
Q_2=-deriv(V,q2);
fprintf('Q_1 = %s  \n',char(Q_1))
fprintf('Q_2 = %s  \n',char(Q_2))
fprintf('\n')


% first Lagrange's equation of motion
Lagrange1 = LHS1-Q1;
% second Lagrange's equation of motion
Lagrange2 = LHS2-Q2;

fprintf('e.o.m:  \n')
pretty(simple(Lagrange1)); 
fprintf(' = 0 \n\n\n')
pretty(simple(Lagrange2)); 
fprintf(' = 0 \n')

data = {L, m, g};
datn = {1 , 1, 9.81};

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

fid = fopen('eomE7_4.m','w+'); 
fprintf(fid,'function dx=eomE7_4(t,x)\n');
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

t0 = 0;  tf = 15; time = [0 tf];

x0 = [pi/18 0 pi/6 0]; 

[ts,xs] = ode45(@eomE7_4, time, x0);

x1 = xs(:,1); 
x2 = xs(:,2);
x3 = xs(:,3); 
x4 = xs(:,4);  

figure(1)
subplot(2,1,1),plot(ts,x1*180/pi,'r'),
xlabel('t (s)'),ylabel('q1 (deg)'),grid,
subplot(2,1,2),plot(ts,x3*180/pi,'b'),
xlabel('t (s)'),ylabel('q2 (deg)'),grid

figure(2)
subplot(2,1,1),plot(x1,x2,'r'),
xlabel('q1'),ylabel('dq1'),grid,
subplot(2,1,2),plot(x3,x4,'b'),
xlabel('q2'),ylabel('dq2'),grid

[tsn,xs] = ode45(@eomE7_4,0:1:5,x0);

fprintf('Results \n\n')
fprintf ...
('    t(s)     q1(rad)  dq1(rad/s) q2(rad)  dq2(rad/s) \n') 
[tsn,xs]

% end of program

% sd={L,m,g,q1,diff(q1,t),q2,diff(q2,t)};
% nd={1,1,9.81,x1,x2,x3,x4};
% Tn=subs(T,sd,nd);
% Vn=double(subs(V,sd,nd));
% 
% figure(3)
% plot(ts,Tn+Vn,'r'),xlabel('t(s)'),ylabel('H(J)'),grid




% Results 
% 
%     t(s)     q1(rad)  dq1(rad/s) q2(rad)  dq2(rad/s) 
% 
% ans =
% 
%          0    0.1745         0    0.5236         0
%     1.0000   -0.2766    0.0877   -0.2067   -1.2662
%     2.0000    0.1416    1.1217    0.2134   -0.4471
%     3.0000    0.0193   -0.2199   -0.1697   -1.8570
%     4.0000    0.0036    0.6631   -0.2805    0.8579
%     5.0000    0.2376   -0.8388    0.1264    0.1029
