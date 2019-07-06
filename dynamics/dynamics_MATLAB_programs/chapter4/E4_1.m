% example 4.1

clear all; clc; close all

syms m g alpha v0 t

N = m*g*cos(alpha);

sx=dsolve('D2x - g*sin(alpha) = 0');
sy=dsolve('D2y = 0');

% dsolve - symbolic solution of ODE
% by default, the independent variable is 't'
% the letter 'D' denotes differentiation d/dt
% 'D' followed by a digit denotes repeated 
%     differentiation
% 'D2' is d^2/dt^2

fprintf('equations of motion: \n')
fprintf('x = %s \n',char(sx))
fprintf('y = %s \n',char(sy))

fprintf ...
('equations of motion with initial conditions\n')
x=dsolve('D2x-g*sin(alpha)=0','Dx(0)=0','x(0)=0');
y=dsolve('D2y = 0','Dy(0) = v0','y(0) = 0');
fprintf('x = %s  \n',char(x))
fprintf('y = %s  \n\n',char(y))

% velocity 
v_x=diff(x,t);
v_y=diff(y,t);
v=[v_x v_y 0];
fprintf('velocity v = \n')
pretty(v); fprintf('\n\n')

% magnitude of velocity
magn_v=sqrt(v(1)^2+v(2)^2)+v(3)^2;
fprintf('|v| = \n')
pretty(magn_v); fprintf('\n\n')

% acceleration
a_x=diff(x,t,2);
a_y=diff(y,t,2);
a_z=0;

a=[a_x a_y a_z];
fprintf('a = \n')
pretty(a); fprintf('\n\n');

% magnitude of acceleration
magn_a=sqrt(a(1)^2+a(2)^2+a(3)^2);
fprintf('|a| = \n')
pretty(simplify(magn_a)); fprintf('\n\n')


% numerical application
gn = 9.81; % m/s^2
alphan = pi/6;
v0n= 4; % m/s

axis manual
axis equal
hold on
grid on

a = 10;
axis([-a a -a a -a a])

x_O=0; y_O=0; z_O=0;
x_B=10; y_B=-10; z_B=0;
x_C=-10; y_C=-10; z_C=0;
x_E=-10; y_E=10; z_E=0;
x_F=10; y_F=10; z_F=0;

t1=text(x_O+0.2, y_O, z_O,' O');
t3=text(x_B, y_B, z_B,' B');
t4=text(x_C, y_C, z_C,' C');

line([x_C,x_B],[y_C,y_B],[z_C,z_B]);
line([x_B,x_F],[y_B,y_F],[z_B,z_F]);
line([x_F,x_E],[y_F,y_E],[z_F,z_E]);
line([x_E,x_C],[y_E,y_C],[z_E,z_C]);

view(10,40); 
light('Position',[1 3 2]); 

% trajectory of the particle
for tn = 0 : 0.005 : 1.5  
   slist={g,alpha,v0,t};
   nlist={gn,alphan,v0n,tn};
   xn = subs(x,slist,nlist);
   yn = subs(y,slist,nlist);
   hm=plot(xn,yn,'k.','Color','red');
   ht=plot(xn,yn);
   pause(0.001)
   delete(hm);
end

% end of program