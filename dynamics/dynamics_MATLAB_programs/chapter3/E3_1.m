% example 3.1
clear all; clc; close all 

syms R omega t;
% parametric equation of motion
theta=omega*t;
x=R*cos(theta);
y=R*sin(theta);

% velocity
v_x=diff(x,t);
v_y=diff(y,t);
v=[v_x v_y];
% magnitude of velocity
magn_v=sqrt(v(1)^2+v(2)^2);

% acceleration
a_x=diff(x,t,2);
a_y=diff(y,t,2);
a=[a_x a_y];
% magnitude of acceleration
magn_a=sqrt(a(1)^2+a(2)^2);

% print symbolic results
fprintf('velocity components v_x and v_y \n')
pretty(v); fprintf('\n\n')
fprintf('magnitude of the velocity \n')
pretty(simplify(magn_v)) 
fprintf('\n\n')

fprintf('acceleration components a_x and a_y')
fprintf('\n')
pretty(a); fprintf('\n\n')
fprintf('magnitude of the acceleration ')
fprintf('\n')
pretty(simplify(magn_a))
fprintf('\n\n')

% numerical application
omegan=1.5;
Rn= 5.0;
slist={R,omega,t};

%set up the axis
axis manual
axis equal
hold on
grid on
axis([-8 8 -8 8])
title('Trajectory, velocity, and acceleration')

%define a scale factor
scalefactor = 3;

for tn = 0 : 0.01 : 3*pi/2 
   nlist = {Rn, omegan, tn};
   xn = subs(x, slist, nlist);
   yn = subs(y, slist, nlist);
   vxn = subs(v_x, slist, nlist)/scalefactor;
   vyn = subs(v_y, slist, nlist)/scalefactor;
   axn = subs(a_x, slist, nlist)/scalefactor;
   ayn = subs(a_y, slist, nlist)/scalefactor;
 
   % position represented by a red dot
   % updated at each tm value of the loop
   hm = plot(xn, yn, 'k.','Color', 'red');
   % trajectory represented by small blue dots
   ht = plot(xn, yn);
    
   % velocity represented by a black vector
   pv=...
    quiver(xn,yn,vxn,vyn,'Color','k','LineWidth',1);
   % acceleration represented by a red vector
   pa=...
    quiver(xn,yn,axn,ayn,'Color','r','LineWidth',1);
   % delay the simulation
   pause(0.001)
   % particle(red dot) deleted 
   % from its old location 
   delete(hm);
   % velocity vector deleted       
   % from its old location    
   delete(pv);
   % acceleration vector deleted 
   % from its old location 
   delete(pa);
end

% for the last location tn = 3*pi/2 
hl = plot(xn,yn,'k.','Color','r');
pause(0.5)
quiver(xn,yn,vxn,vyn,'color','k','LineWidth',1.3);
pause(0.5)
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