% example 6.7
% R-RTR dynamics
% Newton-Euler eom

clear all
close all
clc

L1 = 0.100; % m
L2 = 0.470; % m
L3 = 0.050; % m
AC = 0.300; % m
h =  0.010; % m
d =  0.005; % m
h3 = 0.050; % m
d3 = 0.020; % m

syms t
theta = sym('theta(t)');
omega = diff(theta,t);
alpha = diff(omega,t);

xB = L1*cos(theta);
yB = L1*sin(theta);
rB = [xB yB 0];
vB = diff(rB,t);
aB = diff(vB,t);

xC = AC;
yC = 0;
rC = [xC yC 0];

theta2 = atan((yB-yC)/(xB-xC));
omega2 = diff(theta2,t);
alpha2 = diff(omega2,t);

omega3 = omega2;
alpha3 = alpha2;

alphav = [0 0 alpha];
alpha2v = [0 0 alpha2];
alpha3v = [0 0 alpha3];

xD = xB+L2*cos(theta2);
yD = yB+L2*sin(theta2);

rD = [xD yD 0];
vD = diff(rD,t);
aD = diff(vD,t);

rC1 = rB/2;
rC2 = (rB+rD)/2;
rC3 = rC;

vC1 = diff(rC1,t);
aC1 = diff(vC1,t);
vC2 = diff(rC2,t);
aC2 = diff(vC2,t);
aC3 = diff(rC3,t,2);

% force analysis

g = 9.807;  % m/s^2
rho = 7850; % kg/m^3

m1 = rho*h*d*L1;
m2 = rho*h*d*L2;

m3e = rho*L3*h3*d3;
m3i = rho*L3*h*d;
m3 = m3e - m3i;

IC1 = m1*(L1^2+h^2)/12;
IC2 = m2*(L2^2+h^2)/12;
IC3 = m3e*(L3^2+h3^2)/12-m3i*(L3^2+h^2)/12;
IA = IC1+m1*(L1/2)^2;

G1 = [0 -m1*g 0];
G2 = [0 -m2*g 0];
G3 = [0 -m3*g 0];

syms F23m M23z F03x F03y
F03 = [F03x F03y 0];
F23x = F23m*sin(theta2);
F23y = F23m*cos(theta2);
F23 = [F23x F23y 0];
M23 = [0 0 M23z];
% link 3
eqn3F = -m3*aC3+F03+G3+F23;
eqn3M = -IC3*alpha3v+M23;
eqn3x = eqn3F(1);
eqn3y = eqn3F(2);
eqn3z = eqn3M(3);

syms F12x F12y
F12 = [F12x F12y 0];
% link 2
eqn2F = -m2*aC2+F12-F23+G2;
eqn2M = -IC2*alpha2v+cross(rB-rC2,F12)+...
        cross(rC-rC2,-F23)-M23;
eqn2x = eqn2F(1);
eqn2y = eqn2F(2);
eqn2z = eqn2M(3);

sol32 = solve...
(eqn2x,eqn2y,eqn2z,eqn3z,...
 'F12x','F12y','F23m','M23z');

F12xs = sol32.F12x;
F12ys = sol32.F12y;
F21 = [-F12xs -F12ys 0];

M0 = 1;
w0 = 4;
M = M0*(1-omega/w0);
Mm = [0 0 M];

% link 1
eqn1M = -IA*alphav+cross(rB,F21)+cross(rC1,G1)+Mm;

subout = {diff(theta,t,2), diff(theta,t), theta};
subin = {'ddtheta','x(2)','x(1)'};

eqn1 = subs(eqn1M(3),subout,subin);
NE = solve(eqn1,'ddtheta');
NEeom = char(NE);

% ODE
fid = fopen('eomE6_7_NE.m','w+');
fprintf(fid,'function out = eomE6_7_NE(t,x)\n\n');
fprintf(fid,'out = zeros(2,1);\n');
fprintf(fid,'out(1) = x(2);\n');
fprintf(fid,'out(2) = %s;\n',NEeom);
fclose(fid); 
cd(pwd);

IC = [pi/6 0];
tf = 5;
time = [0 tf];

% integrate
option0 = odeset('RelTol',1e-3,'MaxStep',1e-3);
[t data] = ode45(@eomE6_7_NE,time,IC,option0);

thetan = data(:,1);
omegan = data(:,2);

figure(1)
subplot(2,1,1)
plot(t,thetan)
xlabel('t [s]')
ylabel('theta [rad]')
grid
subplot(2,1,2)
plot(t,omegan);
xlabel('t [s]')
ylabel('omega [rad/s]')
grid

c1 = cos(thetan);
s1 = sin(thetan);

xB = L1*c1;
yB = L1*s1;
rB = [xB yB zeros(length(xB),1)];

% Define phi
theta2n = atan(yB./(xB-xC*ones(length(xB),1)));

% Define position D
xD = xB + L2*cos(theta2n);
yD = yB + L2*sin(theta2n);
rD = [xD yD zeros(length(xB),1)];

figure(2)
hold on
for i = 1:10:length(xB(:,1))
    
    clf
    hold on
    plot([0 xB(i)],[0 yB(i)],'b','LineWidth',5)
    plot([xB(i) xD(i)],[yB(i) yD(i)],'r','LineWidth',5)
    plot(xC,yC,'sk','MarkerSize',10,'LineWidth',10)
    axis([-0.25 0.75 -0.5 0.5])
    grid on
    pause(1/60);

end


% Mmn = M0*(1-omegan/w0);
% figure
% plot(omegan,Mmn)
% xlabel('omega [rad/s]')
% ylabel('M [N m]')

[ts, xs] = ode45(@eomE6_7_NE,0:1:tf,IC);

fprintf('Results \n\n')
fprintf('    t(s)      q(rad)    dq/dt (rad/s) \n'); 
[ts,xs]

% end of program