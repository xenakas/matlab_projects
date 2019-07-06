% example 6.8
% 3D system
% Newton-Euler e.o.m 
clear all; clc; close all

syms t alpha L  m1 m2 g 

q1 = sym('q1(t)');
q2 = sym('q2(t)');

s1 = sin(q1);
c1 = cos(q1);

sa = sin(alpha);
ca = cos(alpha);

% transformation matrix from RF2 to RF1
R21 = [ca sa 0; -sa ca 0; 0 0 1];

% angular velocity of link 1 in RF0 expressed 
% in terms of RF1 {i1,j1,k1}
w1_1 = [diff(q1,t); 0; 0];
fprintf('w1_1 = \n');pretty(w1_1);fprintf('\n')

% angular velocity of link 1 in RF0 expressed 
% in terms of RF2 {i2,j2,k2}
w1_2  = R21*w1_1;
fprintf('w1_2 = \n');pretty(w1_2);fprintf('\n')

% angular acceleration of link 1 in RF0 expressed 
% in terms of RF1{i1,j1,k1}
alpha1_1 = diff(w1_1,t);
fprintf('alpha1_1 = \n');pretty(alpha1_1)
fprintf('\n')

% angular acceleration of link 1 in RF0 expressed 
% in terms of RF2{i2,j2,k2}
alpha1_2 = diff(w1_2,t); 
fprintf('alpha1_2 = \n');pretty(alpha1_2)
fprintf('\n')

% position vector of mass center C1 of link 1 
% in RF0 expressed in terms of RF2 {i2,j2,k2}
rC1_2 = [L; 0; 0];

% position vector of mass center C1 of link 1 
% in RF0 expressed in terms of RF1 {i1,j1,k1}
rC1_1 = R21.'*rC1_2;
fprintf('rC1_1 = \n');pretty(rC1_1)
fprintf('\n')

% linear velocity of mass center C1 of link 1 
% in RF0 expressed in terms of RF1 {i1,j1,k1}
vC1_1 = diff(rC1_1,t) + cross(w1_1, rC1_1);
fprintf('vC1_1 = \n');pretty(vC1_1)
fprintf('\n')

% linear velocity of mass center C1 of link 1 
% in RF0 expressed in terms of RF2 {i2,j2,k2}
vC1_2 = diff(rC1_2,t) + cross(w1_2, rC1_2);
% vC1_2 = simple(R21.'*vC1_1)
fprintf('vC1_2 = \n');pretty(vC1_2);fprintf('\n')

% linear acceleration of mass center C1 of link 1 
% in RF0 expressed in terms of RF1 {i1,j1,k1}
aC1_1 = simple(diff(vC1_1,t)+cross(w1_1,vC1_1));
fprintf('aC1_1 = \n');pretty(aC1_1)
fprintf('\n')

% linear acceleration of mass center C1 of link 1 
% in RF0 expressed in terms of RF2 {i2,j2,k2}
aC1_2 = simple(diff(vC1_2,t)+cross(w1_2,vC1_2));
% aC1_2 = simple(R21.'*aC1_1)
fprintf('aC1_2 = \n');pretty(aC1_2)
fprintf('\n')

% position vector of P 
% in RF0 expressed in terms of RF2 {i2,j2,k2}
rP_2 = [q2; 0;0];
fprintf('rP_2 = \n');pretty(rP_2)
fprintf('\n')

% position vector of mass center C1 of link 1 
% in RF0 expressed in terms of RF1 {i1,j1,k1}
rP_1 = R21.'*rP_2;
fprintf('rP_1 = \n');pretty(rP_1)
fprintf('\n')

% linear velocity of P 
% in RF0 expressed in terms of RF2 {i2,j2,k2}
vP_2 = diff(rP_2,t) + cross(w1_2, rP_2);
fprintf('vP_2 = \n'); pretty(vP_2)
fprintf('\n')

% linear acceleration of P 
% in RF0 expressed in terms of RF2 {i2,j2,k2}
aP_2 = diff(vP_2,t) + cross(w1_2, vP_2);
aP_2 = simple(aP_2);
fprintf('aPx_2 = \n'); pretty(aP_2(1))
fprintf('\n')
fprintf('aPy_2 = \n'); pretty(aP_2(2))
fprintf('\n')
fprintf('aPz_2 = \n'); pretty(aP_2(3))
fprintf('\n')

% gravitational force that acts on link 1 at C1
% RF0 expressed in terms of RF1 {i1,j1,k1}
G1 = [-m1*g; 0; 0] 	

% gravitational force that acts P
% in RF0 expressed in terms of RF2 {i2,j2,k2}
G2 = R21*[-m2*g; 0; 0];
fprintf('G2 = \n'); pretty(G2)
fprintf('\n')

% inertia matrix for link 1 associated with RF2
IC1_2 = [0 0 0; 0 m1*(2*L)^2/12 0; 0 0 m1*(2*L)^2/12];

% F01 reaction force of 0 that acts on link 1 
% at O expressed in terms of RF1 {i1,j1,k1}
syms F01x F01y F01z
F01 = [F01x; F01y; F01z];

% F12 reaction force of 1 that acts on particle 
% at P expressed in terms of RF2 {i2,j2,k2}
syms F12y F12z
F12 = [0; F12y; F12z];

% F21 reaction force of P that acts on link 1 
% at P expressed in terms of RF1 {i1,j1,k1}
F21 = -R21.'*F12;

% reaction moment of 0 that acts on link 1  
% expressed in terms of RF1 {i1,j1,k1}
syms  M01y M01z
M01 = [0; M01y; M01z];

% external torque on link 1
% expressed in terms of RF1 {i1,j1,k1}
syms  T01x 
T01 = [T01x; 0; 0];

% Newton's e.o.m. for link 1: 
% m1 aC1 = F01+G1-F12 =>
eqF1 = F01+F21-m1*aC1_1;
eqF1x = eqF1(1);
eqF1y = eqF1(2);
eqF1z = eqF1(3);

fprintf('F01+G1-F12-m1 aC1 = 0 =>\n\n')
pretty(eqF1x); fprintf('\n')
pretty(eqF1y); fprintf('\n')
pretty(eqF1z); fprintf('\n')
fprintf('\n')

% Newton's e.o.m. for particle P: 
% m2 aP = G2+F12 =>
eqF2 = G2+F12-m2*aP_2;
eqF2x = eqF2(1);
eqF2y = eqF2(2);
eqF2z = eqF2(3);

fprintf('G2+F12-m2*aP = 0 =>\n\n')
pretty(eqF2x); fprintf('\n')
pretty(eqF2y); fprintf('\n')
pretty(eqF2z); fprintf('\n')
fprintf('\n')

% sum of moments about C1 for link 1
M_C1 = cross(-rC1_1,F01) + cross(rP_1-rC1_1,F21) + ...
       M01 + T01;

Min_C1 = -IC1_2*alpha1_2-cross(w1_2,IC1_2*w1_2);
  
% Euler (rotational) e.o.m. for link 1 w.r.t. C1
eqM1 = M_C1+Min_C1;
eqM1x = eqM1(1); 
eqM1y = eqM1(2); 
eqM1z = eqM1(3); 

fprintf('M_C1+Min_C1 = 0 =>\n\n')
pretty(eqM1x); fprintf('\n')
pretty(eqM1y); fprintf('\n')
pretty(eqM1z); fprintf('\n')
fprintf('\n')

ql = {diff(q1,t,2), diff(q2,t,2), ...
   diff(q1,t), diff(q2,t), q1, q2};    
qf = {'ddq1', 'ddq2', ...
   'x(2)', 'x(4)', 'x(1)', 'x(3)'};

sys=...
{eqF1x,eqF1y,eqF1z,eqF2x,eqF2y,eqF2z,eqM1x,eqM1y,eqM1z};     

syseq = subs(sys, ql, qf);

unkn={F01x,F01y,F01z,F12y,F12z,M01y,M01z,'ddq1','ddq2'};

sol=solve(syseq,...
    'F01x,F01y,F01z,F12y,F12z,M01y,M01z,ddq1,ddq2');

ddq1s = simple(sol.ddq1);
ddq2s = simple(sol.ddq2);

% e.o.m
fprintf('d x(1)/dt = x(2)\n')
fprintf('d x(2)/dt = '); pretty(ddq1s);
fprintf('\n')
fprintf('d x(3)/dt = x(4)\n')
fprintf('d x(3)/dt = '); pretty(ddq2s);
fprintf('\n')

data = {alpha, L, m1, m2, g, T01x};
datn = {pi/4, 2, 10, 0.001, 9.81, 300*sin(3*t)};

NE1n = subs(ddq1s, data, datn);
NE2n = subs(ddq2s, data, datn);

% system of ODE 
dx2dt = char(NE1n); 
dx4dt = char(NE2n);

fid = fopen('eomE6_8.m','w+'); 
fprintf(fid,'function dx = eomE6_8(t,x)\n');
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

t0 = 0; tf = 1; time = [0 tf];
x0 = [0 0 1.5 0]; 

 option0 = odeset...
('RelTol',1e-3,'MaxStep',1e-3,'Events',@eventE6_8);

[t,xs,te,ye] = ode45(@eomE6_8, time, x0, option0); 

x1 = xs(:,1); % q1
x2 = xs(:,2); % dq1
x3 = xs(:,3); % q2
x4 = xs(:,4); % dq2  

% plot the results
alpha = pi/4;
xP1 = cos(alpha)*x3;
yP1 = sin(alpha)*x3;
zP1 = 0;

xP0 = x3*cos(alpha);
yP0 = sin(alpha)*x3.*cos(x1);
zP0 = sin(alpha)*x3.*sin(x1);

for kk=1:30:length(xs)    
       
gd = 1.5;
gdm = .25;
axis([-gdm gd -gdm gd -gdm gd]);

set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|-.|--|:')
  
grid on
xlabel('x_0 (m)') 
ylabel('y_0 (m)')
zlabel('z_0 (m)')
title('Position of P')
  
plot3(xP0(kk),zP0(kk),yP0(kk),'r.');
hold on

% Cartesian axes
% quiver3(0,0,0,gd,0,0,1,'Color','k','LineWidth',1)
% text('Interpreter','latex','String','  $x$',...
%     'Position',[gd,0,0],'FontSize',14)
% quiver3(0,0,0,0,gd,0,1,'Color','k','LineWidth',1)
% text('Interpreter','latex','String','  $y$',...
%     'Position',[0,gd,0],'FontSize',14)
% quiver3(0,0,0,0,0,gd,1,'Color','k','LineWidth',1)
% text('Interpreter','latex','String','  $z$',...
%     'Position',[0,0,gd],'FontSize',14)


pause(1/100)

end  

% end of program

