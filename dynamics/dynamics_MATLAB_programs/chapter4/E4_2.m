% example 4.2

clear all; clc; close 

syms m_p v_p m_w v_w x_p x_w L

v_w=solve(m_p*v_p-m_w*v_w,v_w);
fprintf('velocity of the wagon \n')
fprintf('v_w = %s  \n',char(v_w))
fprintf('\n')

m=m_p+m_w;
fprintf('m = %s  \n',char(m))
fprintf('\n')

x_c=(m_p*x_p-m_w*x_w)/m;
fprintf('x_c = %s  \n',char(x_c))
fprintf('\n')

sol=solve(x_c,x_p+x_w-L/2,x_p,x_w);
fprintf('x_p = %s  \n',char(sol.x_p))
fprintf('x_w = %s  \n',char(sol.x_w))
fprintf('\n')

% end of program
