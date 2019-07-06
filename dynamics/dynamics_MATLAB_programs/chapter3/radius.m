clear all; clc; close all

syms t rho

x = sym('x(t)');
y = sym('y(t)');

v = sqrt(diff(x,t)^2+diff(y,t)^2);

a1 = diff(x,t,2)^2+diff(y,t,2)^2;

a2 = diff(v,t)^2+v^4/rho^2;

rhos = simple(solve(a1-a2,rho));

pretty(rhos)



%%

slist = {diff('x(t)',t,2),diff('y(t)',t,2),...
         diff('x(t)',t),diff('y(t)',t)};
nlist = {'ddx', 'ddy', 'dx', 'dy'};

rho1 = subs(rhos,slist,nlist)

rho2 = simple(rho1)

pretty(rho2)

   