% example 3.3
clear all; clc; close all

syms R omega t r real

theta=omega*t;
x=r*(theta-sin(theta));
y=r*(1-cos(theta));

v_x=diff(x,t);
v_y=diff(y,t);
v=[v_x v_y];
magn_v=simplify(sqrt(v(1)^2+v(2)^2));

fprintf('velocity \n')
pretty(v); fprintf('\n\n')
fprintf('magnitude of velocity \n')
pretty(magn_v); fprintf('\n\n')

% position vector of M
rM=[x y];
% position vector of C'
rCp=[r*theta 0];
% position vector C'M
rCpM=simplify(rM-rCp);
rCpMv=simplify(dot(rCpM,v));
fprintf('rCpM.v = %s \n',char(rCpMv))
% rC'M.v = 0 => v perpendicular to rC'M
fprintf('\n\n')

% position vector of C''
rCpp=[r*theta 2*r];
% position vector C''M
rCppM=simplify(rM-rCpp);
rCppM=simplify(dot(rCpM,rCppM));
fprintf('rCpM.rCppM = %s \n',char(rCppM))
% rC'M.rC''M=0=>C'M perpendicular to C''M
fprintf('\n\n')

magn_CpM=sqrt(rCpM(1)^2+rCpM(2)^2);
om=simplify(magn_v/magn_CpM);
fprintf('|v|/|rCpM| = %s \n',char(om))
fprintf('\n\n')

a_x=diff(x,t,2);
a_y=diff(y,t,2);
a=[a_x a_y];
magn_a=sqrt(a(1)^2+a(2)^2);

fprintf('acceleration \n')
pretty(a); fprintf('\n\n')
fprintf('magnitude of acceleration \n')
pretty(simplify(magn_a)); fprintf('\n\n')

% position vector of C
rC=[r*theta r];
rMC=simplify(rC-rM);
a==omega^2*rMC

% simulation of the particle
omegan=2.5;
rn= 0.5;
scale_factor=2.5;
axis manual
axis equal
hold on
grid on
axis([-1 8 -1 2])

for tn = 0 : 0.01 : pi 
   slist={r,omega,t};
   nlist={rn,omegan,tn};
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
   
   [legend_h,object_h,plot_h,text_strings]=...
   legend('Trajectory','Velocity','Acceleration');
   
   set(plot_h(1),'color','b');
   set(object_h(1),'color','b');
   set(plot_h(2),'color','k');
   set(object_h(2),'color','k');
   set(plot_h(3),'color','r');
   set(object_h(3),'color','r');
   
% end of program