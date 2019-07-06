% example 4.5

clear all; clc; close 

syms m g vB k x_0 v r h real;

% kinetic energy T_A of the particle at A 
vA = 0;
T_A = 1/2*m*vA^2;
% kinetic energy T_B of the particle at B
T_B=1/2*m*vB^2;
% potential energy V_A of the particle at A
V_A=1/2*k*x_0^2;
% potential energy V_B of the particle at B
V_B=m*g*(h+r);

vBsol=solve(m*g-m*v^2/r,v);
% velocity of the particle at B 
fprintf('vB =  %s \n',char(vBsol(1)))
fprintf('\n')

T_B=subs(T_B,vB,vBsol(1));
fprintf('T_B =  %s \n',char(T_B))
fprintf('\n')

k=solve(T_A+V_A-(T_B+V_B),k);
fprintf('T_A +V_A = T_B+V_B => \n\n')
fprintf('k =  %s \n',char(vBsol(1)))
fprintf('\n')

ls = {m, h, r, x_0, g};
ln = {1, 0.4, 0.2, 0.08, 9.8};

kn = subs(k, ls, ln);
fprintf('k = %6.3f (N/m)\n', kn)

% end of program