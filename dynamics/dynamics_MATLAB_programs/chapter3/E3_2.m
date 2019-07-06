% example 3.2
clear all; clc; close 

syms R omega t l_1 l_2;

theta=omega*t;
x=l_1*sin(theta);
y=l_2*cos(theta);

v_x=diff(x,t);
v_y=diff(y,t);
v=[v_x v_y];
magn_v=sqrt(v(1)^2+v(2)^2);

fprintf('velocity components v_x and v_y')
fprintf('\n');
pretty(v); fprintf('\n\n')
fprintf('magnitude of velocity')
fprintf('\n');
pretty(magn_v); fprintf('\n\n')

a_x=diff(x,t,2);
a_y=diff(y,t,2);
a=[a_x a_y];
magn_a=sqrt(a(1)^2+a(2)^2);

fprintf('acceleration components a_x and a_y')
fprintf('\n')
pretty(a); fprintf('\n\n')
fprintf('magnitude of acceleration')
fprintf('\n');
pretty(simplify(magn_a)); fprintf('\n\n');

% numerical application
omegan=1.5;
l_1n= 5.0;
l_2n= 2.5;
scale_factor=3;
axis manual
axis equal
hold on
grid on
axis([-6 6 -6 6])
title('Trajectory, velocity, and acceleration')
slist={l_1,l_2,omega,t};
for tn = 0 : 0.01 : 3*pi/2 
   nlist={l_1n,l_2n,omegan,tn};
   xn = subs(x,slist,nlist);
   yn = subs(y,slist,nlist);
   hm=plot(xn,yn,'k.','Color','red');
   ht=plot(xn,yn);
   pause(0.001)
   delete(hm);
end

   hm=plot(xn,yn,'k.','Color','r');
   pause(1.0)
   vxn = subs(v_x,slist,nlist)/scale_factor;
   vyn = subs(v_y,slist,nlist)/scale_factor;
quiver(xn,yn,vxn,vyn,'color','k','LineWidth',1.3);
   pause(1.0)
   axn = subs(a_x,slist,nlist)/scale_factor;
   ayn = subs(a_y,slist,nlist)/scale_factor;
quiver(xn,yn,axn,ayn,'color','r','LineWidth',1.3);
      
[legend_h, object_h, plot_h,text_strings]=...
legend('Trajectory','Velocity','Acceleration');
   set(plot_h(1),'color','b');
   set(object_h(1),'color','b');
   set(plot_h(2),'color','k');
   set(object_h(2),'color','k');
   set(plot_h(3),'color','r');
   set(object_h(3),'color','r');

% end of program