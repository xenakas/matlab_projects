% example 6.4(a)
% dynamics of a double pendulum
% using an inline function

clear all; clc; close all

L1 = 1; L2 = 1; % m
m1 = 1; m2 = 1; % kg
g = 9.8;        % m/s^2

t = sym('t','real');
q1 = sym('q1(t)');
q2 = sym('q2(t)');

xA = L1*cos(q1);
yA = L1*sin(q1);
rA = [xA yA 0];
rC1 = rA/2;
vC1 = diff(rC1,t); 
aC1 = diff(vC1,t);

xB = xA + L2*cos(q2);
yB = yA + L2*sin(q2);
rB = [xB yB 0];
rC2 = (rA + rB)/2;
vC2 = diff(rC2,t);
aC2 = diff(vC2,t);

omega1 = [0 0 diff(q1,t)];
alpha1 = diff(omega1,t);
omega2 = [0 0 diff(q2,t)];
alpha2 = diff(omega2,t);

G1 = [0 -m1*g 0];
G2 = [0 -m2*g 0];

IC1 = m1*L1^2/12;
IO = IC1 + m1*(L1/2)^2;
IC2 = m2*L2^2/12;

% LINK 2
% Sum F for link 2:
% -m2 aC2 + G2 + (-F21) = 0 => 
F21 = -m2*aC2 + G2;

% LINK 1
% Sum M for 1 wrt O:"
% -IO alpha1 + OA x F21 + OC1 x G1 = 0
EqO =...
-IO*alpha1+cross(rA, F21)+cross(rC1,G1);

% LINK 2
% Sum M for 2 wrt C2:
% -IC2 alpha2 + C2A x (-F21) = 0
Eq2 =...
-IC2*alpha2+cross(rA-rC2,-F21);

% list for the symbolical variables 
slist={diff('q1(t)',t,2),diff('q2(t)',t,2),...
 diff('q1(t)',t),diff('q2(t)',t),'q1(t)','q2(t)'};
nlist={'ddq1','ddq2','x(2)','x(4)','x(1)','x(3)'};
% diff('q1(t)',t,2) will be replaced by 'ddq1'
% diff('q2(t)',t,2) will be replaced by 'ddq2'
% diff('q1(t)',t) will be replaced by 'x(2)'
% diff('q2(t)',t) will be replaced by 'x(4)'
% 'q1(t)' will be replaced by 'x(1)'
% 'q2(t)' will be replaced by 'x(3)'

eq1 = subs(EqO(3),slist,nlist);
eq2 = subs(Eq2(3),slist,nlist);

sol = solve(eq1,eq2,'ddq1, ddq2');

dx1 = sym('x(2)');
dx2 = sol.ddq1;
dx3 = sym('x(4)');
dx4 = sol.ddq2;

dx1dt = char(dx1); 
dx2dt = char(dx2);
dx3dt = char(dx3); 
dx4dt = char(dx4);

% digits(3)
fprintf('equations of motion\n') 
fprintf('d[x(1)]/dt=%s\n',char(vpa(dx1dt))) 
fprintf('d[x(2)]/dt=%s\n',char(vpa(dx2dt))) 
fprintf('d[x(3)]/dt=%s\n',char(vpa(dx3dt))) 
fprintf('d[x(4)]/dt=%s\n',char(vpa(dx4dt))) 

g=inline(sprintf('[%s;%s;%s;%s]',...
    dx1dt,dx2dt,dx3dt,dx4dt),'t','x');

t0 = 0;  tf = 10; time = [0 tf];

% initial conditions
x0 = [-pi/4; 0; -pi/4; 0]; 

[t,xs] = ode45(g, time, x0);

x1 = xs(:,1); 
x2 = xs(:,2);
x3 = xs(:,3); 
x4 = xs(:,4);  

subplot(2,1,1),plot(t,x1*180/pi,'r'),...
xlabel('t (s)'),ylabel('q1 (deg)'),grid,...
subplot(2,1,2),plot(t,x3*180/pi,'b'),...
xlabel('t (s)'),ylabel('q2 (deg)'),grid

figure(1)
subplot(2,1,1),plot(t,x1*180/pi,'r'),...
xlabel('t (s)'),ylabel('q1 (deg)'),grid,...
subplot(2,1,2),plot(t,x3*180/pi,'b'),...
xlabel('t (s)'),ylabel('q2 (deg)'),grid

figure(2)
subplot(2,1,1),
plot(x1,x2,'r'),
xlabel('q1'),ylabel('d(q1)/dt'),grid,...
subplot(2,1,2),
plot(x3,x4,'b'),
xlabel('q2'),ylabel('d(q2)/dt'),grid

[ts,xs] = ode45(g,0:1:5,x0);
[ts,xs]
% end of program



% ans =
% 
%          0   -0.7854         0   -0.7854         0
%     1.0000   -2.1833   -1.5597   -2.3810    0.0933
%     2.0000   -1.3296    1.7166   -0.9298    2.0942
%     3.0000   -1.5470   -0.9623   -1.4853   -4.2822
%     4.0000   -1.8816    1.7998   -2.2992    1.5263
%     5.0000   -0.8946   -1.3517   -0.7751    0.2739
