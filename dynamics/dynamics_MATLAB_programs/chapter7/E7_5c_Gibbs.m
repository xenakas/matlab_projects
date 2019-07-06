% example 7.5(c)
% RRT robot arm
% Gibbs-Appell Direct Dynamics

clear; clc; close;

syms t L1 L2 L m1 m2 m3 g real

q1 = sym('q1(t)');
q2 = sym('q2(t)');
q3 = sym('q3(t)');

c1 = cos(q1);
s1 = sin(q1);
c2 = cos(q2);
s2 = sin(q2);

% transformation matrix from RF1 to RF0
R10 = [[1 0 0]; [0 c1 s1]; [0 -s1 c1]];

% transformation matrix from RF2 to RF1
R21 = [[c2 0 -s2]; [0 1 0]; [s2 0 c2]];

% angular velocity of link 1 in RF0 
% expressed in terms of RF1{i1,j1,k1}
w10 = [diff(q1,t) 0 0 ]

% angular velocity of link 2 in RF0 
% expressed in terms of RF1{i1,j1,k1}
w201 = [diff(q1,t) diff(q2,t) 0];

% angular velocity of link 2 in RF0 
% expressed in terms of RF2{i2,j2,k2}
w20 = w201 * transpose(R21) 

% angular acceleration of link 1 in RF0 
% expressed in terms of RF1{i1,j1,k1}
alpha10 = diff(w10,t);

% angular acceleration of link 2 in RF0 
% expressed in terms of RF2{i2,j2,k2}
alpha20 = diff(w20,t); 

% position vector of mass center C1 of link 1 
% in RF0 expressed in terms of RF1{i1,j1,k1}
rC1 = [L1 0 0];

% linear velocity of mass center C1 of link 1 
% in RF0 expressed in terms of RF1{i1,j1,k1}
vC1 = diff(rC1,t) + cross(w10, rC1)

% position vector of mass center C2 of link 2 
% in RF0 expressed in terms of RF1{i1,j1,k1}
rC2 = [L2 0 0];

% linear velocity of mass center C2 of link 2 
% in RF0 expressed in terms of RF1 {i1,j1,k1}
vC2 = simple(diff(rC2,t) + cross(w10,rC2))

% position vector of mass center C3 of link 3 
% in RF0 expressed in terms of RF2{i2,j2,k2}
rC3 = rC2*R21.' + [0 0 q3]

% linear velocity of mass center C3 of link 3 
% in RF0 expressed in terms of RF2{i2,j2,k2}
vC3 = simple(diff(rC3,t) + cross(w20,rC3))

% linear velocity of  C32 of link 2 in RF0
% expressed in terms of RF2{i2,j2,k2} 
% C32 of link 2 is superposed with C3 of link 3
vC32 = simple(vC2 + cross(w20,[0 0 q3]))

% another way of computing vC3 is: 
% vC3p= vC32+diff([0 0 sym('q3(t)')],t);
% vC3-vC3p

% linear accelerations
aC1 = simple(diff(vC1,t)+cross(w10,vC1))
aC2 = simple(diff(vC2,t)+cross(w20,vC2))
aC3 = simple(diff(vC3,t)+cross(w20,vC3));

% kinetic energy
% inertia dyadic

syms I1x I1y I1z I2x I2y I2z real

% inertia matrix associated with 
% central inertia dyadic for link 1 
% expressed in terms of RF1{i1,j1,k1}
I1 = [I1x 0 0; 0 I1y 0; 0 0 I1z];

% inertia matrix associated with 
% central inertia dyadic for link 2
% expressed in terms of RF2{i2,j2,k2}
I2 = [I2x 0 0; 0 I2y 0; 0 0 I2z];

% inertia matrix associated with 
% central inertia dyadic for link 3
% expressed in terms of RF2{i2,j2,k2}
I3 = [m3*L^2/12 0 0; 
      0 m3*L^2/12 0; 
      0 0 0];

% Energy of acceleration
S1 = (1/2)*m1*aC1*aC1.' + (1/2)*alpha10*I1*alpha10.';
S2 = (1/2)*m2*aC2*aC2.' + (1/2)*alpha20*I2*alpha20.';
S3 = (1/2)*m3*aC3*aC3.' + (1/2)*alpha20*I3*alpha20.';
S = expand(S1 + S2 + S3); 

Sddq1 = deriv(S, diff(q1,'t',2));
Sddq2 = deriv(S, diff(q2,'t',2));
Sddq3 = deriv(S, diff(q3,'t',2));

syms T01x T01y T01z T12x T12y T12z F23x F23y F23z real

% contact torque of 0 that acts on link 1 
% in RF0 expressed in terms of RF1{i1,j1,k1}
T01 = [T01x T01y T01z];

% contact torque of link 1 that acts on link 2 
% in RF0 expressed in terms of RF2{i2,j2,k2}
T12 = [T12x T12y T12z];

% contact force of link 2 that acts on link 3 at C3 
% in RF0 expressed in terms of RF2{i2,j2,k2}
F23 = [F23x F23y F23z];

% gravitational force that acts on link 1 at C1
% RF0 expressed in terms of RF1{i1,j1,k1}
G1 = [-m1*g 0 0]	

% gravitational force that acts on link 2 at C2
% in RF0 expressed in terms of RF2{i1,j1,k1}
G2 = [-m2*g 0 0]

% gravitational force that acts on link 3 at C3
% in RF0 expressed in terms of RF2{i2,j2,k2}
G3 = [-m3*g 0 0]*transpose(R21)

% partial velocities 
w1_1 = deriv(w10, diff(q1,t)); 
w2_1 = deriv(w20, diff(q1,t)); 
w1_2 = deriv(w10, diff(q2,t)); 
w2_2 = deriv(w20, diff(q2,t)); 
w1_3 = deriv(w10, diff(q3,t)); 
w2_3 = deriv(w20, diff(q3,t)); 

vC1_1 = deriv(vC1, diff(q1,t)); 
vC2_1 = deriv(vC2, diff(q1,t));
vC1_2 = deriv(vC1, diff(q2,t)); 
vC2_2 = deriv(vC2, diff(q2,t));
vC1_3 = deriv(vC1, diff(q3,t)); 
vC2_3 = deriv(vC2, diff(q3,t));

vC32_1 = deriv(vC32, diff(q1,t)); 
vC3_1 = deriv(vC3, diff(q1,t));
vC32_2 = deriv(vC32, diff(q2,t)); 
vC3_2 = deriv(vC3, diff(q2,t));
vC32_3 = deriv(vC32, diff(q3,t)); 
vC3_3 = deriv(vC3, diff(q3,t));

% generalized active forces

% generalized active force Q1
Q1=w1_1*T01.'+vC1_1*G1.'+w1_1*(R21.')*(-T12.')+...
   w2_1*T12.'+vC2_1*G2.'+vC32_1*(-F23.')+ ...
   vC3_1*F23.'+vC3_1*G3.';

% generalized active force Q2
Q2=w1_2*T01.'+vC1_2*G1.'+w1_2*(R21.')*(-T12.')+...
   w2_2*T12.'+vC2_2*G2.'+vC32_2*(-F23.')+...
   vC3_2*F23.'+vC3_2*G3.';

% generalized active force Q3
Q3=w1_3*T01.'+vC1_3*G1.'+w1_3*(R21.')*(-T12.')+...
   w2_3*T12.'+vC2_3*G2.'+vC32_3*(-F23.')+...
   vC3_3*F23.'+vC3_3*G3.';
 
fprintf('Q1 = %s\n',char(simple(Q1)))
fprintf('Q2 = %s\n',char(simple(Q2))) 
fprintf('Q3 = %s\n',char(simple(Q3)))

% Q1 = T01x
% Q2 = T12y - g*m3*cos(q2(t))*q3(t)
% Q3 = F23z - g*m3*sin(q2(t))
 
% eom
G_A1 = Sddq1-Q1;
G_A2 = Sddq2-Q2;
G_A3 = Sddq3-Q3;

% control torques and control force
q1f=pi/3; q2f=pi/3; q3f=0.3;
b01=450; g01=300;
b12=200; g12=300;
b23=150; g23=50;

T01xc = -b01*diff(q1,t)-g01*(q1-q1f);
T12yc = -b12*diff(q2,t)-g12*(q2-q2f)+g*m3*c2*q3;
F23zc = -b23*diff(q3,t)-g23*(q3-q3f)+g*m3*s2;

tor  = {T01x, T12y, F23z};
torf = {T01xc,T12yc,F23zc};

G_A1 = subs(G_A1, tor, torf);
G_A2 = subs(G_A2, tor, torf);
G_A3 = subs(G_A3, tor, torf);

data = ...
{L1,L2,L,I1x,I2x,I2y,I2z,m1,m2,m3,g};
datn = ...
{0.4,0.4,0.5,5,4,1,4,90,60,40,9.81};

G_A1 = subs(G_A1, data, datn);
G_A2 = subs(G_A2, data, datn);
G_A3 = subs(G_A3, data, datn);

ql = {diff(q1,t,2), diff(q2,t,2), diff(q3,t,2), ...
      diff(q1,t), diff(q2,t), diff(q3,t), q1, q2, q3};    
qf = {'ddq1', 'ddq2',  'ddq3',...
      'x(2)', 'x(4)', 'x(6)', 'x(1)', 'x(3)', 'x(5)'};

% ql                    qf
%----------------------------
% diff('q1(t)',t,2) -> 'ddq1'
% diff('q2(t)',t,2) -> 'ddq2'
% diff('q3(t)',t,2) -> 'ddq3'
%   diff('q1(t)',t) -> 'x(2)'
%   diff('q2(t)',t) -> 'x(4)'
%   diff('q3(t)',t) -> 'x(6)'
%           'q1(t)' -> 'x(1)'
%           'q2(t)' -> 'x(3)' 
%           'q3(t)' -> 'x(5)' 

G_A1 = subs(G_A1, ql, qf);
G_A2 = subs(G_A2, ql, qf);
G_A3 = subs(G_A3, ql, qf);

% solve e.o.m. for ddq1, ddq2, ddq3
sol = solve(G_A1, G_A2, G_A3,'ddq1, ddq2, ddq3');
G_A1 = sol.ddq1;
G_A2 = sol.ddq2;
G_A3 = sol.ddq3;

% system of ODE 
dx2dt = char(G_A1); 
dx4dt = char(G_A2);
dx6dt = char(G_A3);

fid = fopen('eomE7_5c.m','w+'); 
fprintf(fid,'function dx = eomE7_5c(t,x)\n');
fprintf(fid,'dx = zeros(6,1);\n');
fprintf(fid,'dx(1) = x(2);\n');
fprintf(fid,'dx(2) = ');
fprintf(fid,dx2dt);
fprintf(fid,';\n');
fprintf(fid,'dx(3) = x(4);\n');
fprintf(fid,'dx(4) = ');
fprintf(fid,dx4dt);
fprintf(fid,';\n');
fprintf(fid,'dx(5) = x(6);\n');
fprintf(fid,'dx(6) = ');
fprintf(fid,dx6dt);
fprintf(fid,';');
fclose(fid); 
cd(pwd);

t0 = 0;  tf = 15; time = [0 tf];

x0 = [pi/18 0 pi/6 0 0.25 0]; 

[t,xs] = ode45(@eomE7_5c, time, x0);

x1 = xs(:,1); 
x2 = xs(:,2);
x3 = xs(:,3); 
x4 = xs(:,4);  
x5 = xs(:,5); 
x6 = xs(:,6);  

subplot(3,1,1),plot(t,x1*180/pi,'r'),...
xlabel('t (s)'),ylabel('q1 (deg)'),grid,...
subplot(3,1,2),plot(t,x3*180/pi,'b'),...
xlabel('t (s)'),ylabel('q2 (deg)'),grid,...
subplot(3,1,3),plot(t,x5,'g'),...
xlabel('t (s)'),ylabel('q3 (m)'),grid

[ts,xs] = ode45(@eomE7_5c,0:1:5,x0);

fprintf('\n')
fprintf('Results \n'); fprintf('\n');
fprintf...
('t(s) q1 dq1 q2 dq2 q3 dq3 \n'); 
[ts,xs]

% end of program

% Results 
% 
% t(s) q1 dq1 q2 dq2 q3 dq3 
% 
% ans =
% 
%          0    0.1745         0    0.5236         0    0.2500         0
%     1.0000    0.5968    0.3054    0.9316    0.1792    0.2773    0.0196
%     2.0000    0.8184    0.1547    1.0227    0.0380    0.2883    0.0059
%     3.0000    0.9310    0.0787    1.0420    0.0080    0.2925    0.0031
%     4.0000    0.9881    0.0400    1.0461    0.0017    0.2949    0.0019
%     5.0000    1.0172    0.0203    1.0470    0.0004    0.2965    0.0013
