% example 4.3 cartesian

clear all; clc; close all

syms t r alpha m g dtheta

theta = sym('theta(t)');

ur = [ cos(theta), sin(theta), 0];
ut = [-sin(theta), cos(theta), 0];
k = [0, 0, 1];
c = cos(alpha); 
s = sin(alpha);

rP = r*s*ur + r*c*k;
vP = diff(rP,t);

G = [0 0 m*g];
MO = cross(rP, G);
fprintf('MO = r x G = \n')
pretty(MO); fprintf('\n')

HO = cross(rP, m*vP);
fprintf('HO = r x (m v) = \n')
pretty(HO); fprintf('\n')

dHO = diff(HO, t);
fprintf('d(HO)/dt - MO =>\n')
eq = simplify(dHO-MO);
pretty(eq(1)); fprintf('\n')
pretty(eq(2)); fprintf('\n')
pretty(eq(3)); fprintf('\n\n')

ql = {diff(theta,t,2), diff(theta,t)};    
qf = {0, 'dtheta'};

eqf = subs(eq(2), ql, qf);
sol = solve(eqf, 'dtheta');

fprintf('for d2(theta)/dt2 = 0 => d(theta)/dt =\n')
pretty(sol(1))

% end of program