% example 4.10

clear all; clc; close 

syms t m g var omega Ry Rz

x = sym('x(t)');
r = [x 0 0];
omega_v=omega*[cos(var) sin(var) 0];
v_rel=diff(r,t);
a_rel=diff(r,t,2);
a_t=cross(omega_v,cross(omega_v,r));
a_c=2*cross(omega_v,v_rel);
a = a_rel+a_t+a_c;

fprintf('ax = %s\n',char(a(1)))
fprintf('ay = %s\n',char(a(2)))
fprintf('az = %s\n',char(a(3)))
fprintf('\n')

G=-m*g*[cos(var) sin(var) 0];
Rl=[0 Ry Rz];

% eom
eq = m*a-(G+Rl);

fprintf('x =>0 %s=0\n',char(eq(1)))
fprintf('y =>0 %s=0\n',char(eq(2)))
fprintf('z =>0 %s=0\n',char(eq(3)))

% end of program