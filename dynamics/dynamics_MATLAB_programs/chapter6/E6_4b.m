% example 6.4(b)
% dynamics of a double pendulum
% the program uses the function: RR(t,x) 
% the function is defined in the program RR.m
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

F21 = -m2*aC2 + G2;
EqO =...
-IO*alpha1+cross(rA, F21)+cross(rC1,G1);
Eq2 =...
-IC2*alpha2+cross(rA-rC2,-F21);

slist={diff('q1(t)',t,2),diff('q2(t)',t,2),...
 diff('q1(t)',t),diff('q2(t)',t),'q1(t)','q2(t)'};
nlist={'ddq1','ddq2','x(2)','x(4)','x(1)','x(3)'};

eq1 = subs(EqO(3),slist,nlist);
eq2 = subs(Eq2(3),slist,nlist);
sol = solve(eq1,eq2,'ddq1, ddq2');
dx2 = sol.ddq1;
dx4 = sol.ddq2;
dx2dt = char(dx2);
dx4dt = char(dx4);

% opens the file 'RR.m' in the mode specified by 'w+' 
% (create for read and write)
fid = fopen('RR.m','w+');
fprintf(fid,'function dx = RR(t,x)\n');
fprintf(fid,'dx = zeros(4,1);\n');
fprintf(fid,'dx(1) = x(2);\n');
fprintf(fid,'dx(2) = ');
fprintf(fid,dx2dt);
fprintf(fid,';\n');
fprintf(fid,'dx(3) = x(4);\n');
fprintf(fid,'dx(4) = ');
fprintf(fid,dx4dt);
fprintf(fid,';');
% closes the file associated with file identifier fid
fclose(fid);
cd(pwd);
% cd changes current working directory
% pwd displays the current working directory

% fid = fopen(FIL,PERM) opens the file FILE in the
% mode specified by PERM. PERM can be: 
% 'w'     writes (creates if necessary)
% 'w+'    truncates or creates for read and write

t0 = 0;  tf = 10; time = [0 tf];

x0 = [-pi/4 0 -pi/4 0]; % initial conditions

[t,xs] = ode45(@RR, time, x0);

x1 = xs(:,1); 
x3 = xs(:,3); 

subplot(2,1,1),plot(t,x1*180/pi,'r'),...
xlabel('t (s)'),ylabel('q1 (deg)'),grid,...
subplot(2,1,2),plot(t,x3*180/pi,'b'),...
xlabel('t (s)'),ylabel('q2 (deg)'),grid

[ts,xs] = ode45(@RR,0:1:5,x0);
[ts,xs]

% end of program
