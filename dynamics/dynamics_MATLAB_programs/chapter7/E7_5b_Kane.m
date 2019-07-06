% example 7.5(b)
% RRT robot arm
% Kane's dynamical equations

clear; clc; close;

syms t L1 L2 L m1 m2 m3 g real
% generalized coordinates q1, q2, q3
q1 = sym('q1(t)'); 
q2 = sym('q2(t)'); 
q3 = sym('q3(t)');
% generalized speeds u1, u2, u3
u1 = sym('u1(t)'); 
u2 = sym('u2(t)'); 
u3 = sym('u3(t)');
% expressing q1', q2', q3' in terms of 
% generalized speeds u1, u2, u3
dq1 = u1;
dq2 = u2;
dq3 = u3;

qt = {diff(q1,t), diff(q2,t), diff(q3,t)};
qu = {dq1, dq2, dq3};
 
c1 = cos(q1); s1 = sin(q1); 
c2 = cos(q2); s2 = sin(q2);

R10 = [[1 0 0]; [0 c1 s1]; [0 -s1 c1]];
R21 = [[c2 0 -s2]; [0 1 0]; [s2 0 c2]];

w10  = [dq1, 0, 0 ];
w201 = [dq1, dq2, 0];
w20  = w201 * transpose(R21);
alpha10 = diff(w10,t);
alpha20 = subs(diff(w20,t), qt, qu);


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
vC3 = subs(diff(rC3, t), qt, qu) + cross(w20,rC3)

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
aC3 = subs(diff(vC3,t), qt, qu) + cross(w20,vC3)

% the velocities and accelerations are functions of
% q1, q2, q3, u1, u2, u3, and u1', u2', u3'

% partial velocities
w1_1 = deriv(w10, u1); w2_1 = deriv(w20, u1); 
w1_2 = deriv(w10, u2); w2_2 = deriv(w20, u2); 
w1_3 = deriv(w10, u3); w2_3 = deriv(w20, u3);

vC1_1 = deriv(vC1, u1); vC2_1 = deriv(vC2, u1);
vC1_2 = deriv(vC1, u2); vC2_2 = deriv(vC2, u2);
vC1_3 = deriv(vC1, u3); vC2_3 = deriv(vC2, u3);

vC32_1 = deriv(vC32, u1); vC3_1 = deriv(vC3, u1);
vC32_2 = deriv(vC32, u2); vC3_2 = deriv(vC3, u2);
vC32_3 = deriv(vC32, u3); vC3_3 = deriv(vC3, u3);

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

% Kane's dynamical equations

% inertia forces

% inertia force for link 1 
% expressed in terms of RF1{i1,j1,k1}
Fin1= -m1*aC1;
% inertia force for link 2 
% expressed in terms of RF2{i2,j2,k2}
Fin2= -m2*aC2;
% inertia force for link 3 
% expressed in terms of RF2{i2,j2,k2}
Fin3= -m3*aC3;

% inertia moments

% inertia moment for link 1 
% expressed in terms of RF1{i1,j1,k1}
Min1 = -alpha10*I1-cross(w10,w10*I1);
% inertia moment for link 2 
% expressed in terms of RF2{i2,j2,k2}
Min2 = -alpha20*I2-cross(w20,w20*I2);
% inertia moment for link 3 
% expressed in terms of RF2{i2,j2,k2}
Min3 = -alpha20*I3-cross(w20,w20*I3);

% generalized inertia forces

% generalized inertia forces corresponding to q1
Kin1 = w1_1*Min1.' + vC1_1*Fin1.' + ...
       w2_1*Min2.' + vC2_1*Fin2.' + ...
       w2_1*Min3.' + vC3_1*Fin3.';

% generalized inertia forces corresponding to q2
Kin2 = w1_2*Min1.' + vC1_2*Fin1.' + ...
       w2_2*Min2.' + vC2_2*Fin2.' + ...
       w2_2*Min3.' + vC3_2*Fin3.';

% generalized inertia forces corresponding to q3
Kin3 = w1_3*Min1.' + vC1_3*Fin1.' + ...
       w2_3*Min2.' + vC2_3*Fin2.' + ...
       w2_3*Min3.' + vC3_3*Fin3.';

% generalized active forces

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


% generalized active forces 
Q1 = w1_1*T01.' + vC1_1*G1.' + ...
     w1_1*transpose(R21)*(-T12.') + ...
     w2_1*T12.' + vC2_1*G2.' + vC32_1*(-F23.') + ...
     vC3_1*F23.' + vC3_1*G3.';

Q2 = w1_2*T01.' + vC1_2*G1.' + ...
     w1_2*transpose(R21)*(-T12.') + ...
     w2_2*T12.' + vC2_2*G2.' + vC32_2*(-F23.') + ...
     vC3_2*F23.' + vC3_2*G3.';

Q3 = w1_3*T01.' + vC1_3*G1.' + ...
     w1_3*transpose(R21)*(-T12.') + ...
     w2_3*T12.' + vC2_3*G2.' + vC32_3*(-F23.') + ...
     vC3_3*F23.' + vC3_3*G3.';

fprintf('Q1 = %s\n',char(simple(Q1)))
fprintf('Q2 = %s\n',char(simple(Q2))) 
fprintf('Q3 = %s\n',char(simple(Q3)))

% Q1 = T01x
% Q2 = T12y - g*m3*cos(q2(t))*q3(t)
% Q3 = F23z - g*m3*sin(q2(t))
 
% Kane's dynamical equations
% first Kane's dynamical equation
Kane1 = Kin1 + Q1;
% second Kane's dynamical equation
Kane2 = Kin2 + Q2;
% third Kane's dynamical equation
Kane3 = Kin3 + Q3;

% control torques and control force
q1f=pi/3; q2f=pi/3; q3f=0.3;
b01=450; g01=300;
b12=200; g12=300;
b23=150; g23=50;

T01xc = -b01*dq1-g01*(q1-q1f);
T12yc = -b12*dq2-g12*(q2-q2f)+g*m3*c2*q3;
F23zc = -b23*dq3-g23*(q3-q3f)+g*m3*s2;

tor  = {T01x, T12y, F23z};
torf = {T01xc,T12yc,F23zc};

Kan1 = subs(Kane1, tor, torf);
Kan2 = subs(Kane2, tor, torf);
Kan3 = subs(Kane3, tor, torf);

data = ...
{L1,L2,L,I1x,I2x,I2y,I2z,m1,m2,m3,g};
datn = ...
{0.4,0.4,0.5,5,4,1,4,90,60,40,9.81};

Ka1 = subs(Kan1, data, datn);
Ka2 = subs(Kan2, data, datn);
Ka3 = subs(Kan3, data, datn);

ql = {diff(u1,t), diff(u2,t), diff(u3,t) ...
      u1, u2, u3, q1, q2, q3};    
qx = {'du1', 'du2',  'du3',...
'x(4)', 'x(5)', 'x(6)', 'x(1)', 'x(2)', 'x(3)'};

% ql                  qx
%----------------------------
% diff('u1(t)',t) -> 'du1'
% diff('u2(t)',t) -> 'du2'
% diff('u3(t)',t) -> 'du3'
%         'u1(t)' -> 'x(4)'
%         'u2(t)' -> 'x(5)' 
%         'u3(t)' -> 'x(6)' 
%         'q1(t)' -> 'x(1)'
%         'q2(t)' -> 'x(2)' 
%         'q3(t)' -> 'x(3)' 

Du1 = subs(Ka1, ql, qx);
Du2 = subs(Ka2, ql, qx);
Du3 = subs(Ka3, ql, qx);

% solve for du1, du2, du3
sol = solve(Du1, Du2, Du3,'du1, du2, du3');
sdu1 = sol.du1;
sdu2 = sol.du2;
sdu3 = sol.du3;

% system of ODE 
dx1 = char('x(4)');
dx2 = char('x(5)');
dx3 = char('x(6)');
dx4 = char(sdu1);
dx5 = char(sdu2);
dx6 = char(sdu3);

fid = fopen('eomE7_5b.m','w+'); 
fprintf(fid,'function dx = eomE7_5b(t,x)\n');
fprintf(fid,'dx = zeros(6,1);\n');
fprintf(fid,'dx(1) = '); fprintf(fid,dx1); 
fprintf(fid,';\n');
fprintf(fid,'dx(2) = '); fprintf(fid,dx2); 
fprintf(fid,';\n');
fprintf(fid,'dx(3) = '); fprintf(fid,dx3); 
fprintf(fid,';\n');
fprintf(fid,'dx(4) = '); fprintf(fid,dx4); 
fprintf(fid,';\n');
fprintf(fid,'dx(5) = '); fprintf(fid,dx5); 
fprintf(fid,';\n');
fprintf(fid,'dx(6) = '); fprintf(fid,dx6); 
fprintf(fid,';  ');
fclose(fid); cd(pwd);

t0 = 0;  tf = 15; time = [0 tf];

x0 = [pi/18 pi/6 0.25 0 0 0]; 

[t,xs] = ode45(@eomE7_5b, time, x0);

x1 = xs(:,1); 
x2 = xs(:,2);
x3 = xs(:,3); 
x4 = xs(:,4);  
x5 = xs(:,5); 
x6 = xs(:,6);  

subplot(3,1,1),plot(t,x1*180/pi,'r'),...
xlabel('t (s)'),ylabel('q1 (deg)'),grid,...
subplot(3,1,2),plot(t,x2*180/pi,'b'),...
xlabel('t (s)'),ylabel('q2 (deg)'),grid,...
subplot(3,1,3),plot(t,x3,'g'),...
xlabel('t (s)'),ylabel('q3 (m)'),grid

[ts,xs] = ode45(@eomE7_5b,0:1:5,x0);

fprintf('\n')
fprintf('Results \n'); fprintf('\n');
fprintf...
('t(s) q1 q2 q3 u1 u2 u3 \n'); 
[ts,xs]

% end of program