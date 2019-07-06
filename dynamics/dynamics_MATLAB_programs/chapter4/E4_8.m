% example 4.8

clear all; clc; close all

syms m k1 k2 k3 x y z x_0 y_0 z_0 t t_new real

dx_gen=dsolve('m*Dv_x=m*k1');
dy_gen=dsolve('m*Dv_y=m*k2');
dz_gen=dsolve('m*Dv_z=m*k3');

fprintf('velocity \n')
fprintf('dx/dt = %s \n',char(dx_gen))
fprintf('dy/dt = %s \n',char(dy_gen))
fprintf('dz/dt = %s \n',char(dz_gen))

x_gen=dsolve('m*D2x=m*k1');
y_gen=dsolve('m*D2y=m*k2');
z_gen=dsolve('m*D2z=m*k3');

fprintf('position \n');
fprintf('x = %s  \n',char(x_gen))
fprintf('y = %s  \n',char(y_gen))
fprintf('z = %s  \n',char(z_gen))
fprintf('\n')

x_new=dsolve('m*D2x=m*k1','x(0)=x_0','Dx(0)=0');
y_new=dsolve('m*D2y=m*k2','y(0)=y_0','Dy(0)=0');
z_new=dsolve('m*D2z=m*k3','z(0)=z_0','Dz(0)=0');
% initial conditions t=0
% x(0)=x_0, y(0)=y_0, z(0)=z_0
% Dx(0)=0, Dy(0)=0, Dz(0)=0
fprintf('parametric solution of e.o.m. \n')
fprintf('x = %s  \n',char(x_new));
fprintf('y = %s  \n',char(y_new));
fprintf('z = %s  \n',char(z_new));
fprintf('\n')

%replace t^2 by t_new 
x_newr=subs(x_new,t^2,t_new);
y_newr=subs(y_new,t^2,t_new);
z_newr=subs(z_new,t^2,t_new);

t_x = solve(x-x_newr,t_new)/2;
t_y = solve(y-y_newr,t_new)/2;
t_z = solve(z-z_newr,t_new)/2;

fprintf('implicit equation of trajectory \n');
fprintf('%s = %s = %s \n',...
    char(t_x),char(t_y),char(t_z));

axis manual
axis equal
hold on
grid on
a=20;
axis([-a a -a a -a a])
view(30,20) 

for t_n = 0 : 0.02 : 3 
   slist={x_0,y_0,z_0,k1,k2,k3,t};
   nlist={3,3,7,4.0,-1.5,-5.0,t_n};
   x_t=subs(x_new,slist,nlist);
   y_t=subs(y_new,slist,nlist);
   z_t=subs(z_new,slist,nlist);
   hm=plot3(x_t,y_t,z_t,'.','Color','red');
   ht=plot3(x_t,y_t,z_t,'.');
   title('trajectory')
   pause(0.001)
   delete(hm);
 
end

% end of program