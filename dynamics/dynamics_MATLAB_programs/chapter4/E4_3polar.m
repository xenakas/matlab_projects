% example 4.3 polar coordinates

clear all; clc; close all

syms t r alpha m g dtheta

theta = sym('theta(t)');

c = cos(alpha); 
s = sin(alpha);

rP = [r*s, 0, r*c];
omega = [0, 0, diff(theta,t)]; 
v = diff(rP,t) + cross(omega, rP);
fprintf('v = \n')
pretty(v); fprintf('\n')

G = [0 0 m*g];
MO = cross(rP, G);
fprintf('MO = r x G = \n')
pretty(MO); fprintf('\n')

HO = cross(rP, m*v);
fprintf('HO = r x (m v) = \n')
pretty(HO); fprintf('\n')

dHO = diff(HO, t)+cross(omega, HO);
fprintf('d(HO)/dt - MO =>\n')
eq = dHO-MO;
pretty(eq(1)); fprintf('\n')
pretty(eq(2)); fprintf('\n')
pretty(eq(3)); fprintf('\n\n')

eqf = subs(eq(2), diff(theta, t), dtheta);
sol = solve(eqf, 'dtheta');

fprintf('for d2(theta)/dt2 = 0 => d(theta)/dt =\n')
pretty(sol(1))

% end of program