function [PP] = Rotate(P,theta,a,b,c)
% function [PP] = Rotate(P,theta,index)
% Rotate a point (x,y,z) at the angle theta ccw
% about x-axis, if index = 1, 
% about y-axis, if index = 2, and 
% about z-axis, if index = 3
% P - coordinate matrix of an old view, 
% PP = R*P - a new coordinate matrix of a new view

[m,n] = size(P); % m = 3 for all examples

%input the theta_x, theta_y and theta_z angle

omega=1;
theta_x=a/sqrt(a^2+b^2+c^2);
theta_y=b/sqrt(a^2+b^2+c^2);
theta_z=c/sqrt(a^2+b^2+c^2);

%calculate the vector omega 
omega_x=omega*theta_x;
omega_y=omega*theta_y;
omega_z=omega*theta_z;

R=[...
omega_x*omega_x+(1-omega_x*omega_x)*cos(theta) ...
omega_x*omega_y*(1-cos(theta))-omega_z*sin(theta) ...
omega_x*omega_z*(1-cos(theta))+omega_y*sin(theta);...
omega_x*omega_y*(1-cos(theta))+omega_z*sin(theta) ...
omega_y*omega_y+(1-omega_y*omega_y)*cos(theta) ...
omega_y*omega_z*(1-cos(theta))-omega_x*sin(theta);...
omega_x*omega_z*(1-cos(theta))-omega_y*sin(theta) ...
omega_y*omega_z*(1-cos(theta))+omega_x*sin(theta) ...
omega_z*omega_z+(1-omega_z*omega_z)*cos(theta);];      
PP = R*P;
