function [v_r,a_r,magn_v,magn_a] = part_der(x,y,z,t) 

% 1st Order Partial Derivative

% velocity components 
v_x=diff(x,t);
v_y=diff(y,t);
v_z=diff(z,t);
v_r=[v_x v_y v_z];
% magnitude of the velocity
magn_v=sqrt(v_r(1)^2+v_r(2)^2+v_r(3)^2);


% 2nd Order Partial Derivative
a_x=diff(x,t,2);
a_y=diff(y,t,2);
a_z=diff(z,t,2);
a_r=[a_x a_y a_z];
% magnitude of the acceleration
magn_a=sqrt(a_r(1)^2+a_r(2)^2+a_r(3)^2);