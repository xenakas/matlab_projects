% example 7.4(b)
% Hamilton e.o.m.
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

% kinetic energy 
T1=IO*omega1*omega1.'/2;
T2=m2*vC2*vC2.'/2+IC2*omega2*omega2.'/2;
T = T1 + T2;
% potential energy V
V=-m1*g*rC1(2)-m2*g*rC2(2);

% Lagrangian L
Lagrangian = expand(T-V);

p1L = deriv(Lagrangian, diff(q1,t));
p2L = deriv(Lagrangian, diff(q2,t));

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
Hdq=p1L*diff(q1,t)+p2L*diff(q2,t)-Lagrangian;
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

fprintf('Hamilton equations \n\n')
fprintf('dq1/dt = %s \n',char(Hdq1))
fprintf('dq2/dt = %s \n',char(Hdq2))
fprintf('\n')
fprintf('dp1/dt = %s \n',char(Hdp1))
fprintf('dp2/dt = %s \n',char(Hdp2))

data = {L, m, g};
datn = {1 , 1, 9.81};

Hdq1n = subs(Hdq1,data,datn);
Hdq2n = subs(Hdq2,data,datn);
Hdp1n = subs(Hdp1,data,datn);
Hdp2n = subs(Hdp2,data,datn);

slist = {q1, q2, p1, p2};
nlist = {'x(1)','x(2)','x(3)','x(4)'};
eom1 = subs(Hdq1n,slist,nlist);
eom2 = subs(Hdq2n,slist,nlist);
eom3 = subs(Hdp1n,slist,nlist);
eom4 = subs(Hdp2n,slist,nlist);

dx1 = char(eom1);
dx2 = char(eom2);
dx3 = char(eom3);
dx4 = char(eom4);

fH=inline(sprintf...
('[%s;%s;%s;%s]',dx1,dx2,dx3,dx4),'t','x');

t0 = 0;  tf = 15; time = [0 tf];

q10  = pi/18;
q20  = pi/6 ;
dq10 = 0;
dq20 = 0;

qlist  = {diff(q1,t),q1,diff(q2,t),q2};
qlist0 = {dq10,q10,dq20,q20};

p10s = subs(p1L,data,datn);
p10  = double(subs(p10s,qlist,qlist0)); 
p20s = subs(p2L,data,datn);
p20  = double(subs(p20s,qlist,qlist0)); 

x0 = [q10; q20; p10; p20];  
[ts,xs] = ode45(fH, time, x0);

x1 = xs(:,1); 
x2 = xs(:,2);
x3 = xs(:,3); 
x4 = xs(:,4);  

figure(1)
subplot(2,1,1),plot(ts,x1*180/pi,'r'),
xlabel('t (s)'),ylabel('q1 (deg)'),grid,
subplot(2,1,2),plot(ts,x2*180/pi,'b'),
xlabel('t (s)'),ylabel('q2 (deg)'),grid

sd={L,m,g,q1,q2,p1,p2};
nd={1,1,9.81,x1,x2,x3,x4}; 
dx1=double(subs(dq1p,sd,nd));
dx2=double(subs(dq2p,sd,nd));

figure(2)
subplot(2,1,1),plot(x1,dx1,'r'),
xlabel('q1'),ylabel('dq1'),grid,
subplot(2,1,2),plot(x2,dx2,'b'),
xlabel('q2'),ylabel('dq2'),grid

% end of program

% [ts,xs] = ode45(fH,0:1:5,x0);
% fprintf('Results \n\n')
% fprintf ...
% (' t q1 q2 p1 p2\n') 
% [ts,xs]
 
% ans =
% 
%          0    0.1745    0.5236         0         0
%     1.0000   -0.2765   -0.2070   -0.5145   -0.3785
%     2.0000    0.1421    0.2125    1.2727    0.4104
%     3.0000    0.0197   -0.1704   -1.2062   -0.7261
%     4.0000    0.0037   -0.2807    1.2938    0.6054
%     5.0000    0.2367    0.1282   -1.0687   -0.3810

% Tn=double(subs(T,sd,nd));
% Vn=double(subs(V,sd,nd));
% figure(3)
% plot(ts,(Tn+Vn),'r'),xlabel('t(s)'),
% ylabel('H(J)'),grid



